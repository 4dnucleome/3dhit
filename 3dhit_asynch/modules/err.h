/*
 * err.h
 *
 * autor:  Lukasz Bieniasz-Krzywiec
 *
 */

#ifndef ERR_H_
#define ERR_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <string.h>

/*
 * Wypisuje informacje o błędnym zakończeniu funkcji systemowej i kończy
 * działanie.
 */
extern void syserr(const char *fmt, ...);

/*
 * Wypisuje informacje o błędzie i kończy działanie.
 */
extern void fatal(const char *fmt, ...);

#endif /* ERR_H_ */
