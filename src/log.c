// Creates and maintains a log file.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>

// Log file:
void createLogFile(){
  FILE * logFile = fopen("log/log.txt","w+"); // Create the log file. Remove if exists.
  if (logFile == NULL){
    printf("ERROR: Cannot create log file.");
    exit(EXIT_FAILURE);
  }
  else{
    fclose(logFile);
  }
}

// Add a new log entry
void addLogEntry(char *str){
  // Logging variables:
  time_t now = time(0);
  struct tm * timenow = localtime(&now);
  char timeInString[50];

  FILE * logFile = fopen("log/log.txt","a+"); // Append to log file.
  strftime(timeInString, sizeof(timeInString), "%Y-%m-%d %H:%M:%S", timenow); // Format time string.

  fprintf(logFile, "%s\t%s\n", timeInString, str); // Write to log file.
  fclose(logFile);
}
