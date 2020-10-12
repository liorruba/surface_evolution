// Creates and maintains a log file.
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include <fstream>
#include <chrono>
#include "../include/log.hpp"

// Log file:
void createLogFile(const char *path){
  std::ofstream logFile;
  logFile.open(path, std::ios_base::app);
  if(!logFile) {
    printf("ERROR: Cannot create log file.");
    exit(EXIT_FAILURE);
  }
  else {
    addLogEntry("Log created successfully.");
  }
}

// Add a new log entry
void addLogEntry(const char *str){
  std::ofstream logFile;
  // Logging variables:
  logFile.open("log/log.txt", std::ios_base::app); // Append to log file.

  // Logging variables:
  time_t now = time(nullptr);
  logFile << std::put_time(localtime(&now), "%F_%T") << "\t" << str << std::endl; // Write to log file.
}
