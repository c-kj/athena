#pragma once

// 用于在多个自定义文件中 include 的声明

#include <map>

// SourceTerm （主要是 SN 和 cooling）的注入时机。按照单个 cycle 内的先后顺序排列，从而可以进行 < 比较
// 使用 .at 方法来获取对应的 SourceTermPosition，如果不存在则会抛出异常
enum class SourceTermPosition {InSourceTerm, AfterSourceTerm, UserWorkInLoop}; 

// 这里也许不能声明为 static?
static std::unordered_map<std::string, SourceTermPosition> source_term_position_map = {
  {"InSourceTerm", SourceTermPosition::InSourceTerm},
  {"AfterSourceTerm", SourceTermPosition::AfterSourceTerm},
  {"UserWorkInLoop", SourceTermPosition::UserWorkInLoop}
};
