// 这个头文件的目的是声明一些可以在多个 pgen 中复用的内容。

#pragma once

#include <unordered_map>

/* -------------------------------------------------------------------------- */
/*               以下内容普遍适用于各种 pgen，可以在多个 pgen 之间复用。              */
/* -------------------------------------------------------------------------- */

// SourceTerm （主要是 SN 和 cooling）的注入时机。按照单个 cycle 内的先后顺序排列，从而可以进行 < 比较
// 有关源项的时机，参看我的笔记： 《黑洞吸积 SN 湍流 project.md》和《黑洞吸积 SN 湍流 TODO 历史.md》
enum class SourceTermPosition {InSourceTerm, AfterSourceTerm, UserWorkInLoop}; 

// 用于处理 input 参数中的字符串。对这个 map 使用 .at 方法来获取对应的 SourceTermPosition，如果不存在则会抛出异常
//? 这里也许不能声明为 static?
static std::unordered_map<std::string, SourceTermPosition> source_term_position_map = {
  {"InSourceTerm", SourceTermPosition::InSourceTerm},
  {"AfterSourceTerm", SourceTermPosition::AfterSourceTerm},
  {"UserWorkInLoop", SourceTermPosition::UserWorkInLoop}
};