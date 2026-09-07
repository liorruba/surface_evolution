#pragma once
#include <string>

// Creates and maintains a log file.
void createLogFile(const char *path);

// Add a new log entry, optionally echoed to the screen.
void addLogEntry(const char *str, bool dispOnScreen = false);
void addLogEntry(const std::string &str, bool dispOnScreen = false);
