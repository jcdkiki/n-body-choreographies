#include "logging.hpp"
#include "simulation.hpp"

const int LOGGED_BODY_INDEX = 0;

std::ofstream explicit_log;
std::ofstream implicit_log;
std::ofstream adams_log;
std::ofstream adaptive_log;

void init_logs() {
    explicit_log.open("explicit.log");
    implicit_log.open("implicit.log");
    adams_log.open("adams.log");
    adaptive_log.open("adaptive.log");
}

void close_logs() {
    explicit_log.close();
    implicit_log.close();
    adams_log.close();
    adaptive_log.close();
}

void log_position(std::ofstream& log, double t) {
    log << t << " " 
        << state[LOGGED_BODY_INDEX].position.x << " "
        << state[LOGGED_BODY_INDEX].position.y << " "
        << state[LOGGED_BODY_INDEX].position.z << "\n";
}
