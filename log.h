#ifndef LOG_H
#define LOG_H

int open_log(const char *fname);
int flog(const char *fmt, ...);

//extern int loglevel;
extern int echolog;
extern FILE *logfp;

#endif
