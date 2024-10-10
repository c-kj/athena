
#include <sstream>  // std::stringstream, std::getline
#include <fstream>  
#include <iomanip>  // std::setprecision
#include <stdexcept>  // std::invalid_argument
#include <algorithm>  // std::sort
#include <random>  // 生成随机数

// #include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"


#include "supernova.hpp"
#include "utils.hpp"
#include "initial_condition.hpp"  // rho_from_n_cgs 用于从初始条件中的 n_cgs 换算为 rho，用于计算 SN Rate
#include "my_outputs.hpp"


/* -------------------------------------------------------------------------- */
/*                               SupernovaEvent                               */
/* -------------------------------------------------------------------------- */


// 构造函数
SupernovaEvent::SupernovaEvent(SupernovaParameters *paras, Real t, Point x, Vector v): 
  paras(paras), time(t), position(x), velocity(v)
{
  auto region_name = paras->region_name;
  if (region_name == "ball") {
    // energy 的 region 是默认的，input file 里的 radius 和 center 都是指 energy 的。mass 的相关参数可以单独指定。
    energy_region = std::unique_ptr<Region>(new Ball(position, paras->radius, paras->ndim));
    // mass 的 region。如果没有单独指定，则使用 SN_i_radius 和 SN_i_center 所指定的。
    // 目前不检查 mass 是否为 0，避免麻烦。内存开销应该很小。
    mass_region = std::unique_ptr<Region>(new Ball(position, paras->mass_radius, paras->ndim));
    
  } else if (region_name == "heart") {
    //FUTURE 把 heart SN 实现
    throw std::invalid_argument("Heart region is not implemented yet.");
  } else {
    throw std::invalid_argument("Unknown region type: " + region_name);
  }

  // energy 是 input file 中的数字，单位为 SN_energy_unit (1e51 ergs)。volume 是 code unit。因此 energy_density 是 code units。mass 同理。
  energy_density = paras->energy * paras->pSNe->SN_energy_unit / energy_region->volume; //* volume 有时会被初始化为 0，因为这种时候 region.contains() 永远不会为 true
  mass_density   = paras->mass   * paras->pSNe->SN_mass_unit   / mass_region->volume;
}


std::string SupernovaEvent::Info() const {
  std::ostringstream oss;
  oss << "time = " << time << ", position = {" << position << "}, velocity = {" << velocity << "}";
  return oss.str();
}



/* -------------------------------------------------------------------------- */
/*                             SupernovaParameters                            */
/* -------------------------------------------------------------------------- */


// SupernovaParameters 的构造函数
SupernovaParameters::SupernovaParameters(Supernovae *pSNe, ParameterInput *pin, const int i): pSNe(pSNe), ndim(pSNe->pmy_mesh->ndim) {
  //TODO 等把 input file 的 block 格式改写后，可以把这里都改为 const 成员变量，用初始化列表初始化。
  const std::string block_name = "supernova";
  const std::string i_str = std::to_string(i);

  // 初始化自身的参数
  energy      = pin->GetOrAddReal(block_name, "SN_" + i_str + "_energy", 0.0);
  mass        = pin->GetOrAddReal(block_name, "SN_" + i_str + "_mass", 0.0);
  region_name = pin->GetOrAddString(block_name, "SN_" + i_str + "_region", "ball"); // region 的种类名，默认是 ball
  radius      = pin->GetOrAddReal(block_name, "SN_" + i_str + "_radius", 0.0);
  mass_radius = pin->GetOrAddReal(block_name, "SN_" + i_str + "_mass_radius", radius);  // mass 注入区域的半径。默认和 energy radius 一致

  // 准备好用于创建所属的 events 的参数
  std::vector<Real> time_list;
  std::vector<Point> center_list;
  //FUTURE 目前对速度的处理只是统一的平动速度。以后可以考虑加入随机速度、Kepler 速度。
  Vector velocity = read_array(pin->GetOrAddString(block_name, "SN_" + i_str + "_velocity", "0,0,0"));


  std::string event_type = pin->GetOrAddString(block_name, "SN_" + i_str + "_event_type", "single_point"); // 默认是 single_point，即只有一个位置
  if (event_type == "random") {  //* 如果是 random，目前按照 tau_SF 以及区域内总质量来计算时间间隔，等间隔地生成 time_list。注意目前 time_list 并不 random！

    std::string allow_region_shape = pin->GetOrAddString(block_name, "SN_" + i_str + "_allow_region_shape", "ball");
    std::unique_ptr<Region> allow_region_outer;
    std::uniform_real_distribution<Real> dis_x, dis_y, dis_z;
    RegionSize mesh_size = pSNe->pmy_mesh->mesh_size;
    Point BH_pos = {0,0,0};

    // 确定 SN 的生成区域，从而计算爆炸数目（根据体积）、生成所需的随机分布
    if (allow_region_shape == "ball") { // 以黑洞为中心的球体区域
      Real allow_region_radius = pin->GetReal(block_name, "SN_" + i_str + "_allow_region_radius");
      allow_region_outer = std::unique_ptr<Region>(new Ball(BH_pos, allow_region_radius, ndim));  // 注意这里的 allow_region_outer 实际上并不是字面意义的 allow_region，还要减去 sink_region 才行
      // 检验 allow_region_radius 是否超出了模拟的 mesh 的范围，如果超出的话，后续体积计算会变得不正确，因此报错退出。
      std::vector<Real> BH_mesh_boundary_distances = { // 这里的写法不太优雅，但没办法，mesh_size 本身太不优雅了。
          std::abs(mesh_size.x1min - BH_pos[0]), std::abs(mesh_size.x1max - BH_pos[0]),
          std::abs(mesh_size.x2min - BH_pos[1]), std::abs(mesh_size.x2max - BH_pos[1]),
          std::abs(mesh_size.x3min - BH_pos[2]), std::abs(mesh_size.x3max - BH_pos[2])};
      auto max_distance = *std::max_element(BH_mesh_boundary_distances.begin(), BH_mesh_boundary_distances.end());
      if (allow_region_radius > max_distance) {
        throw std::invalid_argument("allow_region_radius > xmax is not allowed!");
      }
      //* generate_region_radius 的用途是，固定 random_seed 的情况下，可以手动控制各个 SN 的位置是随着 allow_region_radius 缩放，随着 xmax 缩放，还是不变
      Real generate_region_radius = pin->GetOrAddReal(block_name, "SN_" + i_str + "_generate_region_radius", allow_region_radius);  // 这个参数一般是不需要显式指定的，默认与 allow_region_radius 相同
      dis_x = std::uniform_real_distribution<Real>(- generate_region_radius, generate_region_radius);
      dis_y = std::uniform_real_distribution<Real>(- generate_region_radius, generate_region_radius);
      dis_z = std::uniform_real_distribution<Real>(- generate_region_radius, generate_region_radius);

    } else if (allow_region_shape == "all") { // 如果是 all，则生成区域是整个网格。用于均匀 ISM 的模拟
      allow_region_outer = std::unique_ptr<Region>(new Cuboid({mesh_size.x1min, mesh_size.x2min, mesh_size.x3min}, {mesh_size.x1max, mesh_size.x2max, mesh_size.x3max}));
      dis_x = std::uniform_real_distribution<Real>(mesh_size.x1min, mesh_size.x1max);
      dis_y = std::uniform_real_distribution<Real>(mesh_size.x2min, mesh_size.x2max);
      dis_z = std::uniform_real_distribution<Real>(mesh_size.x3min, mesh_size.x3max);

    } else {
      throw std::invalid_argument("Unknown region type: " + allow_region_shape);
    }

    // 限制不要让 SN 在黑洞内生成
    Real R_in = pin->GetOrAddReal("problem","R_in", 0.0);  //FUTURE 目前直接从 pin 读取，但以后可以改为从统一的某个类指针里读取。
    Ball sink_region(BH_pos, R_in, ndim);  // 目前暂时让黑洞的位置是原点，以后可能考虑变化，那就要用单独的类来管理了。
    // 允许 SN 爆炸的区域的体积
    Real allow_region_volume = allow_region_outer->volume - sink_region.volume;  
    //* 这里假定了 sink_region 完全包含于 allow_region_outer 中，否则体积计算是不正确的。而且这里是理想的球体体积，而非网格离散化后的体积。
    // 如果想要更通用的体积计算，需要实现 Region 的交集运算 & 交集的体积计算，比较困难。一种方式是用 Monte Carlo 方法。

    // 计算 SN 的爆炸时间 time_list
    {
      Real t_start = pin->GetOrAddReal(block_name, "SN_" + i_str + "_t_start", 0.0); 
      Real tlim = pin->GetReal("time", "tlim"); 
      Real t_end = pin->GetOrAddReal(block_name, "SN_" + i_str + "_t_end", tlim);
      Real tau_SF = pin->GetReal(block_name, "SN_" + i_str + "_tau_SF") * pSNe->punit->million_yr_code;  // 从 input file 中读取 tau_SF，单位是 Myr，转换为 code unit
      //TODO 这里目前从 pin 读取的参数，并重新进行换算。以后考虑改为从 initial_condition 指针中读取。但这需要在顶层把指针收集起来。
      Real n_init_cgs = pin->GetReal("initial_condition", "n_init_cgs");
      Real rho_init_code = rho_from_n_cgs(Abundance::mu, n_init_cgs) / pSNe->punit->code_density_cgs;
      
      // 计算 SN 的爆炸速率。以下都是在 code unit 下计算
      Real SF_rate_density = rho_init_code / tau_SF;  // 恒星形成率密度，量纲：质量密度/时间
      Real SN_rate_density = SF_rate_density / (150. * pSNe->punit->solar_mass_code);  // SN 爆炸率密度，量纲：数密度/时间。每 150 个 Msun 形成，对应一个 SN 爆炸，这个数字来自 Stellar IMF 的积分
      Real SN_rate = SN_rate_density * allow_region_volume;  // SN 爆炸率，量纲：个数/时间
      Real SN_time_interval = 1.0 / SN_rate;  // SN 爆炸的平均时间间隔，量纲：时间
      
      // 生成等间隔的 time_list (相当于 numpy.arange)
      for (Real t = t_start; t < t_end; t += SN_time_interval) {
        time_list.push_back(t);
      }
    }

    // 生成 SN 的位置 center_list
    int random_seed = pin->GetOrAddInteger(block_name, "SN_" + i_str + "_random_seed", 1234);
    std::mt19937 gen(random_seed);
    while (center_list.size() < time_list.size()) {
      Point center = {dis_x(gen), dis_y(gen), dis_z(gen)};  //* 注意这里消耗 gen 的顺序，是 x1, y1, z1, x2, y2, z2, ...，而非 x1, x2, ... y1, y2, ... z1, z2, ... 的顺序
      if (allow_region_outer->contains(center) and not sink_region.contains(center)) {  // 拒绝采样，要求生成的点在 allow_region_outer 内，且不在 sink_region 内
        center_list.push_back(center);
      }
    }

  } else if (event_type == "single_point") {  // 单点 SN 爆炸。时间上可以重复多次。
    time_list = read_vector(pin->GetString(block_name, "SN_" + i_str + "_time"));  // 读取 time_list，可能有多个，也可能只有一个。

    std::string center_str = pin->GetOrAddString(block_name, "SN_" + i_str + "_center", "0,0,0"); // 用于后续的 fallback
    Point center = read_array(center_str);
    center_list.assign(time_list.size(), center);  // 把 center_list 初始化为 event_num 个相同的 center

  } else if (event_type == "from_file") { // 从文件中读取
    throw std::invalid_argument("Event type 'from_file' is not implemented yet.");

  } else {
    throw std::invalid_argument("Unknown event type: " + event_type);
  }

  // 创建所属的 SupernovaEvent 
  std::size_t event_num = time_list.size();
  for (int i = 0; i < event_num; ++i) {  //TODO 这里叫 i 跟函数参数重了，会不会有问题？
    event_list.reserve(event_num);
    event_list.emplace_back(new SupernovaEvent(this, time_list[i], center_list[i], velocity));
    // 使用 emplace_back 避免了 push_back 的移动开销。隐含调用 vector<T> 中 T 的构造函数，也更简洁。
  }

}






/* -------------------------------------------------------------------------- */
/*                                 Supernovae                                 */
/* -------------------------------------------------------------------------- */


// Supernovae 的构造函数
Supernovae::Supernovae(Mesh *pmy_mesh, ParameterInput *pin) : 
    pmy_mesh(pmy_mesh), punit(pmy_mesh->punit),
    SN_energy_unit(punit->bethe_code), SN_mass_unit(punit->solar_mass_code)  // 这里暂时没有写开关：什么情况下把 SN 的 input 单位解读为 1e51 erg 和 Msun。若有需求再改
{
  SN_flag = pin->GetOrAddInteger("supernova","SN_flag",0);  //* 目前 SN_flag 设为一个 int，是为了以后可能会有多种 SN 的情况（暂时也没想出来）
  if (SN_flag == 0) return; // 如果不开启 Supernovae 机制，直接跳出
  if (SN_flag < 0) { throw std::invalid_argument("SN_flag < 0 is not allowed!"); }

  if ( !NON_BAROTROPIC_EOS ) {
    throw std::runtime_error("SN explosion is not implemented for barotropic EOS!"); 
  }
  
  InitSupernovaParameters(pin);
  GatherSupernovaEvents();
  // if (Globals::my_rank == 0) { PrintInfo(); }
  if (Globals::my_rank == 0) { 
    std::cout << "SupernovaEvent Info:\n";
    std::cout << "Total number: " << supernova_list.size() << "\n";
    WriteInfoCSV("info/Supernovae.csv");
  }
}


// 在 Supernovae 初始化时调用，从 pin 中读取各类超新星的参数，存入 supernova_paras_list
void Supernovae::InitSupernovaParameters(ParameterInput *pin) {
  std::string block_name, i_str;
  block_name = "supernova";
  for (int i = 1; true; i++) {
    i_str = std::to_string(i);
    if (pin->DoesParameterExist("supernova", "SN_" + i_str + "_time")) { // 目前根据是否有 SN_i_time 来判断是否有这个 SN
      supernova_paras_list.emplace_back(new SupernovaParameters(this, pin, i));
    } else {
      break;
    }
  }
}


void Supernovae::GatherSupernovaEvents() {
  // 把所有 SupernovaParameters 中的 SupernovaEvent 指针放入 supernova_list
  for (const auto& paras : supernova_paras_list) {  // 遍历所有 SupernovaParameters
    for (const auto& SN: paras->event_list) {  // 把每个 paras 对应的 events（的裸指针）都放入 supernova_list
      supernova_list.emplace_back(SN.get());
    }
  }

  // 把 supernova_list 按时间升序排序
  std::sort(supernova_list.begin(), supernova_list.end(),
            [](const SupernovaEvent *a, const SupernovaEvent *b) {
              return a->time < b->time;
            });
}


void Supernovae::PrintInfo() const {
  //TODO 比较临时，有待完善: 每个 paras 输出：能量、tau_SF 等参数
  std::cout << "Supernovae: \n" 
  << "Total: " << supernova_list.size() << "\n";
  for (const auto& SN : supernova_list) {
    std::cout << SN->Info() << "\n";
  }
}


//TODO 可以考虑改进格式，从而在读取时完全不丢失精度
void Supernovae::WriteInfoCSV(const std::string &filename) const {
  ensure_parent_directory_exists(filename);
  std::ofstream file{filename};

  // 一些总体信息？

  // 所有 SupernovaEvent 的信息
  file << "SupernovaEvent Info:\n";
  file << "Total number: " << supernova_list.size() << "\n";
  // 写入 header
  file << "=========================================================\n";
  file << "t, x, y, z, v_x, v_y, v_z, energy_density, number_density\n";

  // 设置浮点数输出格式，确保不丢失精度
  file << std::fixed << std::setprecision(std::numeric_limits<Real>::max_digits10);

  // 写入每个 SupernovaEvent 的信息
  for (const auto& SN : supernova_list) {
    file << SN->time << ", " << SN->position << ", " << SN->velocity << ", " << SN->energy_density << ", " << SN->mass_density << "\n";
  }
}


void Supernovae::GetSupernovaeToInject(const Real time, const Real dt) {
  // 对比 time 和 dt 是否有变化，避免重复计算。在每个 MPI Rank 上只需要计算一次，而无需对每个 MeshBlock 都计算一次。
  static Real last_time = -1.0, last_dt = -1.0;    // static 初始化为 -1，避免一开始就相等
  if (time == last_time && dt == last_dt) return;  // 如果 time 和 dt 没有变化，直接跳出
  last_time = time; last_dt = dt;

  // 双指针法：寻找 start 和 end，这两个 iterator 用于指示 supernova_list 中当前时间步内需要注入的 SN 的范围。
  const Real end_time = time + dt;
  static auto start = supernova_list.begin();  // start 是 static 的，继承上次的位置，避免每次都从头开始找
  if (time < last_time) { start = supernova_list.begin(); }  // 如果时间倒流 (dt < 0) 或者在 InSourceTerm 时机注入且各 stage 的起始时间并非单调不减，则改为从头开始找。不过这种情形目前应该不会出现。
  // 找到新的 start 和 end。如果是左闭右开 ( time <= t and t < time + dt )，则这两个 lambda 判断都是 >=；若是左开右闭，则都是 >。
  //BUG 这里似乎应该左开右闭？SourceTerm 末尾的时间究竟应该算时间步的末尾还是开头？？？ 
  start = std::find_if(start, supernova_list.end(), 
                      [time](const SupernovaEvent *SN) { return SN->time >= time; }
                      );
  auto end = std::find_if(start, supernova_list.end(), 
                      [end_time](const SupernovaEvent *SN) { return SN->time >= end_time; }
                      );  // 找 end 时，仍然从 start 开始找，因为没准儿 dt 会比上次缩短，end 未必总是单增的。

  supernova_to_inject.assign(start, end); // 把找到的 SN 指针放入 supernova_to_inject。assign 会自动清空，不必手动 clear。
}

void Supernovae::SuperNovaeSourceTerm(MeshBlock *pmb, const Real time, const Real dt,
            const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
            const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
            AthenaArray<Real> &cons_scalar) {

  const Real mesh_time = pmb->pmy_mesh->time;
  const Real mesh_dt = pmb->pmy_mesh->dt;

  // 当前调用 SourceTerm 时的注入比例。
  // 对于 AfterSourceTerm 和 UserWorkInLoop，传入的 dt 都是 mesh_dt，从而 inject_ratio = 1.0
  // 对于 InSourceTerm，传入的 dt 是当前 stage 的 dt，inject_ratio = beta_last_stage，根据当前 stage 的比例来调整注入的比例。
  // 注意：传入到这个函数的 dt 不一定是 总 SourceTerm 函数在当前 stage 接收到的 dt。因此这个变量不叫 dt_ratio，以示区分。
  const Real inject_ratio = dt/mesh_dt;

  GetSupernovaeToInject(mesh_time, mesh_dt); // 寻找应在当前步中爆炸的 SN 时，使用 mesh_time 和 mesh_dt，而非 SourceTerm 传入的 time 和 dt
  if ( supernova_to_inject.empty() ) return; // 如果当前时间步没有超新星爆炸，直接跳出
  
  Real injected_energy = 0.0;
  Real injected_mass   = 0.0;

  // 把 SN 的遍历放在外层，方便最内层 SIMD 矢量化，内存访问也更连续。
  for (SupernovaEvent* SN : supernova_to_inject) {
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      Real z = pmb->pcoord->x3v(k);
      for (int j = pmb->js; j <= pmb->je; ++j) {
        Real y = pmb->pcoord->x2v(j);
  #pragma omp simd reduction(+:injected_energy, injected_mass)  // SIMD 矢量化。这里的循环体较为简单，应该可以成功？但不是性能热点，没多大提升。
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real x = pmb->pcoord->x1v(i);
          Real cell_volume = pmb->pcoord->GetCellVolume(i,j,k);

          // 超新星爆炸的能量和质量注入。
          if ( SN->energy_region->contains({x,y,z}) ) {
            Real dE = SN->energy_density * inject_ratio;  // E 是能量密度而非能量！
            cons(IEN,k,j,i) += dE;
            injected_energy += dE * cell_volume;
          }
          if ( SN->mass_region->contains({x,y,z}) ) {
            Real drho = SN->mass_density * inject_ratio;
            Real dE = 0.5 * drho * (SQR(SN->velocity[0]) + SQR(SN->velocity[1]) + SQR(SN->velocity[2]));  // 注入的质量密度对应的动能密度
            cons(IDN,k,j,i) += drho;
            cons(IM1,k,j,i) += drho * SN->velocity[0];
            cons(IM2,k,j,i) += drho * SN->velocity[1];
            cons(IM3,k,j,i) += drho * SN->velocity[2];
            cons(IEN,k,j,i) += dE;
            injected_mass   += drho * cell_volume;
            injected_energy += dE   * cell_volume;
            // 注入 Passive Scalar，正比于质量密度
            if (NSCALARS > 0) {
              cons_scalar(PassiveScalarIndex::SN,k,j,i) += drho;
            }
          }
          // 这里的 debug 信息不再那么有用，将来可以考虑删去
          // if (debug >= DEBUG_Main && pmb->gid == 0 && SN_flag > 0){
          //   if(mesh_time <= SN_time && SN_time < mesh_time + mesh_dt && dt == mesh_dt ) {
          //     printf_and_save_to_stream(debug_stream, 
          //     "DEBUG_Mesh: SN explosion at time = %f, mesh_time = %f, dt = %f, mesh_dt = %f \n", 
          //     time, pmb->pmy_mesh->time, dt, pmb->pmy_mesh->dt);
          //   }
          // } 
        }
      }
    }
  }
  namespace idx = RealUserMeshBlockDataIndex;
  // 记录 SN 注入的能量、质量、个数
  //* 只有当 SN 的 SourceTerm 在 AfterSourceTerm，也即最后一个 stage 中调用时，以下记录才是准确的。否则由于 integrator 的 多个 stage，要偏大一些。
  pmb->ruser_meshblock_data[idx::SN_injected_energy](0) += injected_energy;
  pmb->ruser_meshblock_data[idx::SN_injected_mass  ](0) += injected_mass;
  pmb->ruser_meshblock_data[idx::SN_injected_number](0) += supernova_to_inject.size(); // 在每个 MeshBlock 上都记录 SN 的个数，各个 MeshBlock 都一样，在 hst 函数中用 min 汇总。

  return;
}


#if false
Real Supernovae::SupernovaeTimeStep(MeshBlock *pmb) const {
  const Real mesh_time = pmb->pmy_mesh->time;
  const Real mesh_dt = pmb->pmy_mesh->dt;
  const Real end_time = mesh_time + mesh_dt;
  
  auto it = std::upper_bound(supernova_list.begin(), supernova_list.end(), end_time,
                             [](const SupernovaEvent *SN, const Real t) {
                               return SN->time < t;
                             });
  if ( it != supernova_list.end() ) {
    //TODO 
    //BUG 这里目前没有考虑开始 t=0 时的情况。
    return (*it)->time - end_time;
  }
  return std::numeric_limits<Real>::max();
}
#endif