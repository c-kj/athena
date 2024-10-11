

// C++ headers
#include <sstream>

// 自定义函数所需的标准库导入
#include <vector>
#include <array>
#include <tuple>

#include <fstream>  // std::ofstream
#include <iostream> // std::cout, std::endl
#include <cstdarg>  // va_list, va_start, va_end, vsprintf 等用于处理变长参数列表
#include <sys/stat.h>  // 为了使用 stat 和 mkdir 函数




// Athena++ headers
#include "../../athena.hpp"
#include "../../parameter_input.hpp"

#include "utils.hpp"


// 自定义的用于将 debug 信息同时 print 和输出到文件的函数
// 定义一个函数，它接收一个输出流、一个字符串前缀、一个格式化字符串和一些可变参数
//! 注意！不知道这个函数如果用于 MPI 并行计算时，会不会出问题！可能会抢占式写入？？？最好是只用在本地 debug 时
//! 对于 MeshBlock 级别以上，OpenMP or MPI 时也是危险的，同时往一个流里面写入
void printf_and_save_to_stream(std::ostream &stream, const char *format, ...) {
    va_list args;           // 定义一个 va_list 类型的变量，用于存储可变参数列表
    va_start(args, format); // 使用 va_start 函数初始化 args，使其包含从 format 之后的所有参数

    char buffer[256]; // 定义一个字符数组，用于存储格式化后的字符串

    // 使用 vsnprintf 函数将格式化的字符串写入 buffer，vsnprintf 函数会根据 format 和 args 生成一个字符串，并确保不会超过 buffer 的大小
    vsnprintf(buffer, sizeof(buffer), format, args);

    printf("%s", buffer); // 使用 printf 函数将 buffer 输出到屏幕

    // 检查 stream 是否有效
    if (stream)
    {
        // 如果 stream 有效，将 output 输出到 stream
        stream << buffer;
    }
    else
    {
        // 如果 stream 无效，输出一条错误消息
        printf("stream is not open !!!!! \n");
    }

    // 使用 va_end 函数结束可变参数列表
    va_end(args); 
}


// 用于检查并创建目录的函数
void check_and_create_directory(const std::string& path) {
  struct stat info;  // 存储文件或目录的信息
  if (stat(path.c_str(), &info) != 0) {
    mkdir(path.c_str(), 0755);  // 如果路径不存在，创建它
  } else if (!(info.st_mode & S_IFDIR)) {
    throw std::runtime_error(path + " exists but is not a directory");  // 如果路径存在，但不是一个目录，抛出错误
  }
}

// 检查一个路径的父目录是否存在，如果不存在则创建它。
// 会创建 path 中最后一个 / 之前的全部路径，因此结尾是否有 / 是不同的
void ensure_parent_directory_exists(const std::string& path) {
  std::string partial_path;  // 存储路径的部分内容

  for (char c : path) {
    partial_path += c;  // 将当前字符添加到部分路径中
    if (c == '/') {
      check_and_create_directory(partial_path);  // 检查并创建部分路径
    }
  }
}

// 用于从 input file 中读取自定义的 AMR 参数，以 point_1 = 0, 1, 2.0 这种形式给出点的坐标
std::tuple<std::vector<std::array<Real, 3>>, std::vector<int>> get_AMR_points_and_levels(ParameterInput *pin) {
  std::string block_name, point_name, level_name, num;
  std::array<Real, 3> point;
  std::vector<std::array<Real, 3>> point_list;
  std::vector<int> level_list;
  int level, default_level;

  default_level = pin->GetInteger("mesh", "numlevel") - 1; // numlevel 是从 1 开始的，所以要减 1
  block_name = "AMR/point";                    // 从 input file 中读取的 Block 的名字

  for (int i = 1; true; i++) {                // 让 i 递增，直到找不到 point_i 为止
    point_name = "point_" + std::to_string(i);
    level_name = "level_" + std::to_string(i);
    if (pin->DoesParameterExist(block_name, point_name)) {
      // 把逗号分隔的三个数字读入 point 数组
      point = read_array(pin->GetString(block_name, point_name));
      level = pin->GetOrAddInteger(block_name, level_name, default_level);   // 把 level_i 读入 level
      // 把 point 和 level 添加到 point_list 和 level_list 中
      point_list.push_back(point);
      level_list.push_back(level);

    } else {
      break;
    }
  }
  return std::make_tuple(point_list, level_list);    // 返回一个 tuple
}

// 读取 3 个数组成的列表，用于 Point、Vector 等
std::array<Real, 3> read_array(std::string str, char delimiter) {
  int size = 3;
  std::array<Real, 3> array;
  std::vector<Real> vec = read_vector(str, delimiter);
  if (vec.size() != size) { // 如果解析出的列表长度与 array 所要求的长度不符，抛出异常
    throw std::invalid_argument("The input string should contain exactly 3 numbers separated by " + std::string(1, delimiter));
  }
  for (int j = 0; j < size; j++) {
    array[j] = vec[j];
  }
  return array;
}

//TODO 可以考虑用 src/utils/string_utils.cpp 中的 split 函数重写？也许功能更强一些。
// 从字符串中读取数字列表
// 能解析 "1"、"1,2,3"、"1,2,3,"，但不能 "1,2,3,,"
std::vector<Real> read_vector(std::string str, char delimiter) {
  std::string num;
  std::vector<Real> vec;
  
  std::stringstream ss(str);
  while (std::getline(ss, num, delimiter)) {
    vec.push_back(std::stod(num));
  }
  return vec;
}

// 返回表示时间间隔的字符串
std::string format_duration(const std::chrono::duration<double>& duration) {
  using namespace std::chrono;
  std::ostringstream oss;
  // 自定义 days 类型
  using days = std::chrono::duration<int, std::ratio<86400>>;

  auto d = duration_cast<days>(duration);
  auto h = duration_cast<hours>(duration - d);
  auto m = duration_cast<minutes>(duration - d - h);
  auto s = duration_cast<seconds>(duration - d - h - m);

  if (d.count() > 0) {
    oss << d.count() << "d ";
  }
  if (h.count() > 0 || d.count() > 0) {
    oss << h.count() << "h ";
  }
  if (m.count() > 0 || h.count() > 0 || d.count() > 0) {
    oss << m.count() << "m ";
  }
  oss << s.count() << "s";

  return oss.str();
}