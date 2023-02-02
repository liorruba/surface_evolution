// Creates and maintains a log file.
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include <fstream>
#include <chrono>
#include <iomanip>
#include "../include/log.hpp"

// Log file:
void createLogFile(const char *path){
  std::ofstream logFile;
  logFile.open(path, std::ios_base::app);
  if(!logFile) {
    std::cout << "ERROR: Cannot create log file." << std::endl;
    exit(EXIT_FAILURE);
  }
  else {
    addLogEntry("Log created successfully.", false);
  }
}

// Add a new log entry
void addLogEntry(const char *str, bool dispOnScreen = false){
  std::ofstream logFile;
  // Logging variables:
  logFile.open("log/log.txt", std::ios_base::app); // Append to log file.

  // Logging variables:
  time_t now = time(nullptr);
  logFile << std::put_time(localtime(&now), "%F_%T") << "\t" << str << std::endl; // Write to log file.

  if (dispOnScreen) std::cout << std::put_time(localtime(&now), "%F_%T") << "\t" << str << std::endl;
}
