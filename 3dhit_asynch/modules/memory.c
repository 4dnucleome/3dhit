/*
 * tools.c
 *
 * autor: Łukasz Bieniasz-Krzywiec
 *
 */

#include "../3dhit.h"

/*
 * Alokuje/realokuje pamięć.
 */
void *mem(void *p, size_t size) {
  p = realloc(p, size);
  if (p == NULL) {
    fatal("Cannot alocate %u bytes of memory.\n", size);
  }
  return p;
}
