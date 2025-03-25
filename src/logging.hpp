
#pragma once
#include <fstream>

extern std::ofstream explicit_log;
extern std::ofstream implicit_log;
extern std::ofstream adams_log;
extern std::ofstream adaptive_log;

void init_logs();
void close_logs();
void log_position(double t);
