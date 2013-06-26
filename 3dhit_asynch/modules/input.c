/*
 * input.c
 *
 * autor: Łukasz Bieniasz-Krzywiec
 *
 */

#include "../3dhit.h"

/*
 * Wczytuje białko z pliku.
 */
prot_t *readpdb(char *pdbfn, char ***xpdb) {
  const char num2seq[21] = {
    'G','A','S','C','V','T','I','P','M','D',
    'N','L','K','E','Q','R','H','F','Y','W','X'
  };

  const char num2name[21][4] = {
    "GLY","ALA","SER","CYS","VAL","THR","ILE","PRO",
    "MET","ASP","ASN","LEU","LYS","GLU","GLN","ARG",
    "HIS","PHE","TYR","TRP","GAP"
  };

  int           npdb = 0, i, j, k;
  unsigned int  len = 0, len16 = 0;
  char          line[MAX_LINE_LENGTH + 1];

  line[MAX_LINE_LENGTH] = '\0';

  FILE *pdbfile = fopen(pdbfn, "rt");

  if (pdbfile == NULL) {
    fatal("There is no %s PDB file.\n", pdbfn);
  }

  while (fgets(line, MAX_LINE_LENGTH, pdbfile) != NULL) {
    npdb += 1;
    if (!(strncmp(line + 13, "CA", 2) || strncmp(line, "ATOM", 4))) {
      len += 1;
    }
  }

  fclose(pdbfile);

  len16 = len + 16 - MODULO(len, 4);

  prot_t *prot;
  posix_memalign((void*)(&prot), 128, sizeof(char*) * npdb);

  char **pdb = (char**)mem(NULL, sizeof(char*) * npdb);

  prot->len = len;
  prot->name = pdbfn;
  posix_memalign((void*)(&prot->num), 128, sizeof(short) * len16);
  posix_memalign((void*)(&prot->seq), 128, sizeof(char) * (len16 + 1));
  posix_memalign((void*)(&prot->coor), 128, sizeof(coor_t) * len16);
  posix_memalign((void*)(&prot->c), 128, sizeof(coor_t) * len16);

  pdbfile = fopen(pdbfn, "rt");

  j = k = 0;
  while (fgets(line, MAX_LINE_LENGTH, pdbfile) != NULL) {
    i = strlen(line);
    pdb[j] = (char*)malloc((i + 1) * sizeof(char));
    strncpy(pdb[j], line, i);
    pdb[j][i] = '\0';

    j += 1;

    if (strncmp(line + 13, "CA", 2) || strncmp(line, "ATOM", 4)) {
      continue;
    }

    for (i = 0; i < 20; ++i) {
      if(!strncmp(line + 17, num2name[i], 3)) break;
    }

    prot->seq[k] = num2seq[i];
    prot->num[k] = atoi(line + 22);
    prot->coor[k][0] = atof(line + 30);
    prot->coor[k][1] = atof(line + 38);
    prot->coor[k][2] = atof(line + 46);
    line[13] = 'X';

    k += 1;
  }

  prot->seq[len] = '\0';
  *xpdb = pdb;

  return prot;
}
