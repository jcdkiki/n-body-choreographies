#include "logging.hpp"
#include "simulation.hpp"

const int LOGGED_BODY_INDEX = 0;

std::ofstream logfile;

void init_logs() {
    logfile.open("file.log");
}

void close_logs() {
    logfile.close();
}

void log_position(double t) {
    logfile << t << " " 
        << state[LOGGED_BODY_INDEX].position.x << " "
        << state[LOGGED_BODY_INDEX].position.y << " "
        << state[LOGGED_BODY_INDEX].position.z << "\n";
}
