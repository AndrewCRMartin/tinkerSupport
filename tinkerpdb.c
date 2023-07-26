/*************************************************************************

   Program:    tinkerpdb
   File:       tinkerpdb
   
   Version:    V1.0
   Date:       17.09.15
   Function:   Convert a Tinker .xyz file into a PDB file - the Tinker
               xyzpdb program seems to be unreliable!
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2015
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew.martin@ucl.ac.uk
               andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

   ONLY WORKS WITH THE AMBER99 PARAMETER FILE

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   17.09.15  Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/fsscanf.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF        240
#define MAXATOMTYPES  5000
#define MAXTYPELABEL    32
#define MAXLABEL         8
#define MAXWORDS         8
#define MAXCHAINLABEL    8
#define CNDISTSQ       3.5
#define CADISTSQ      16.0
#define TINKERDATA    "TINKERDATA"

/* Used to store information about fields to look up from the Tinker
   parameter file
*/
typedef struct _lookup
{
   char *input1;
   int  input1Len;
   char *input2;
   int  input2Len;
   int  outField;
   char *output;
   int  atomField;
   BOOL het;
   int  nfields;
}  LOOKUP;


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *paramFile, 
                  char *infile, char *outfile, char ***chains);
BOOL tinker2pdb(FILE *in, FILE *paramFp, char **chains, FILE *out);
void ReadTinkerAtomTypes(FILE *fp, char resnam[MAXATOMTYPES][MAXLABEL], 
                         char atnam[MAXATOMTYPES][MAXLABEL],
                         BOOL isHet[MAXATOMTYPES], int maxTypes);
void Usage(void);
void ExtractTypesFromTinkerAtomRecord(char *buffer, char *resnam, 
                                      char *atnam, BOOL *isHet);
BOOL ConvertTinkerDescriptionToResnamAndAtnam(
   char words[MAXWORDS][MAXTYPELABEL], int nwords,
   char *resnam, char *atnam, BOOL *isHet);
void PopulatePDBRecord(PDB *p, int atnum, REAL x, REAL y, REAL z,
                       char *resnam, char *atnam, BOOL isHet);
PDB *ReadTinkerAsPDB(FILE *in, FILE *paramFp, char *header);
void FixHydrogens(PDB *pdb);
void FixCterOxygens(PDB *pdb);
void FixAtomNames(PDB *pdb);
void FixILECD1(PDB *pdb);
void InsertNumberInAtnam(PDB *q, int hydrogenNumber);
char *GetChainLabel(int ChainNum);
void DoChain(PDB *pdb, char **chains, BOOL BumpChainOnHet);
void doFixAtomName(PDB *start, PDB *stop, char *atom, int nChars);
void RenumberResidues(PDB *pdb);


/************************************************************************/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        paramFile[MAXBUFF],
        **chains = NULL;
   FILE *in  = stdin,
        *out = stdout,
        *pFp = NULL;
   BOOL noEnv = FALSE;
   
    
   if(ParseCmdLine(argc, argv, paramFile, infile, outfile, &chains))
   {
      if((pFp=blOpenFile(paramFile, TINKERDATA, "r", &noEnv))==NULL)
      {
         fprintf(stderr,"Error: Unable to open Tinker parameter \
file: %s\n", paramFile);
         if(noEnv)
         {
            fprintf(stderr,"       Try setting %s environment \
variable.\n", TINKERDATA);
         }
         
         return(1);
      }
      
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if(!tinker2pdb(in, pFp, chains, out))
         {
            fprintf(stderr,"Error: Conversion failed\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"Error: Unable to open input of output file\n");
         return(1);
      }
      
   }
   else
   {
      Usage();
   }
   
   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *paramFile, 
                     char *infile, char *outfile, char ***chains)
   --------------------------------------------------------------
   Input:   int    argc          Argument count
            char   **argv        Argument array
   Output:  char   *paramFile
            char   *infile       Input file (or blank string)
            char   *outfile      Output file (or blank string)
            char   ***chains     Chain labels
   Returns: BOOL                 Success?

   Parse the command line

   17.09.15  Original   By: ACRM  
*/
BOOL ParseCmdLine(int argc, char **argv, char *paramFile, 
                  char *infile, char *outfile, char ***chains)
{
   argc--;
   argv++;

   infile[0]   = outfile[0] = paramFile[0] = '\0';
   
   if(argc < 1)
   {
      return(FALSE);
   }
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         if(argv[0][2]!='\0')
         {
           return(FALSE);
         }
         else
         {            
            switch(argv[0][1])
            {
            case 'h':
               return(FALSE);
               break;
            case 'c':
               if(!(--argc))
                  return(FALSE);

               argv++;

               if((*chains=blSplitStringOnCommas(argv[0], MAXCHAINLABEL))
                  ==NULL)
               {
                  fprintf(stderr,"No memory for storing chain labels: \
%s\n",
                          argv[0]);
                  exit(1);
               }
               
               break;
            default:
               return(FALSE);
               break;
            }
         }         
      }
      else
      {
         /* Check that there are 1, 2 or 3 arguments left               */
         if(argc < 1 || argc > 3)
            return(FALSE);
         
         /* Copy the first to paramFile                                 */
         strcpy(paramFile, argv[0]);
         argc--;
         argv++;

         /* If there's another, copy it to infile                       */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
            argv++;
         }
         
         /* If there's another, copy it to outfile                      */
         if(argc)
         {
            strcpy(outfile, argv[0]);
            argc--;
            argv++;
         }
         
         return(TRUE);
      }
      
      argc--;
      argv++;
   }
   return(TRUE);
}

/************************************************************************/
BOOL tinker2pdb(FILE *in, FILE *paramFp, char **chains, FILE *out)
{
   PDB  *pdb;
   char header[MAXBUFF];
   
   if((pdb=ReadTinkerAsPDB(in, paramFp, header))==NULL)
      return(FALSE);

   /* Apply chain labels                                                */
   DoChain(pdb, chains, FALSE);
   RenumberResidues(pdb);

   WritePDB(out, pdb);

   return(TRUE);
}

/************************************************************************/
void ReadTinkerAtomTypes(FILE *fp, char resnam[MAXATOMTYPES][MAXLABEL], 
                         char atnam[MAXATOMTYPES][MAXLABEL],
                         BOOL isHet[MAXATOMTYPES], int maxTypes)
{
   char buffer[MAXBUFF],
        atomType[MAXTYPELABEL];
   int  atnum;

   /* Clear the output arrays                                           */
   for(atnum=0; atnum<MAXATOMTYPES; atnum++)
   {
      resnam[atnum][0] = '\0';
      atnam[atnum][0]  = '\0';
   }

   /* Read the 'atom' records from the file                             */
   while(fgets(buffer, MAXBUFF, fp))
   {
      if(!strncmp(buffer,"atom   ", 7))
      {
         fsscanf(buffer,"%10x%5d%15x%27s",&atnum, atomType);
         ExtractTypesFromTinkerAtomRecord(atomType, 
                                          resnam[atnum], atnam[atnum],
                                          &(isHet[atnum]));

#ifdef DEBUG
         fprintf(stdout, "%5d \"%-4s\" : \"%-4s\"\n", 
                 atnum, resnam[atnum], atnam[atnum]);
#endif
      }
   }
}

/************************************************************************/
void ExtractTypesFromTinkerAtomRecord(char *buffer, 
                                      char *resnam, 
                                      char *atnam,
                                      BOOL *isHet)
{
   char *ptr,
         words[MAXWORDS][MAXTYPELABEL];
   int   nWords;

   /* Remove the double inverted commas                                 */
   ptr = buffer+1;
   TERMAT(ptr, '"');

   /* Blank the words                                                   */
   for(nWords=0; nWords<MAXWORDS; nWords++)
      words[nWords][0] = '\0';

   /* Split the string into words                                       */
   nWords=0;
   while((ptr=blGetWord(ptr, words[nWords], MAXTYPELABEL))!=NULL)
   {
      if(nWords>=MAXWORDS)
         break;
      nWords++;
   }
   nWords++;

#ifdef DEBUG
   {
      int i;
      for(i=0; i<MAXWORDS; i++)
      {
         fprintf(stdout,"%s :", words[i]);
      }
      fprintf(stdout,"\n");
   }
#endif

   ConvertTinkerDescriptionToResnamAndAtnam(words, nWords, 
                                            resnam, atnam, isHet);
}




/************************************************************************/
BOOL ConvertTinkerDescriptionToResnamAndAtnam(
   char words[MAXWORDS][MAXTYPELABEL], int nwords,
   char *resnam, char *atnam, BOOL *isHet)
{
   static LOOKUP
      lookup[] = {{"Gly",      3, NULL,    0, 0, "GLY ", 1, FALSE, 2},
                  {"Ala",      3, NULL,    0, 0, "ALA ", 1, FALSE, 2},
                  {"Val",      3, NULL,    0, 0, "VAL ", 1, FALSE, 2},
                  {"Leu",      3, NULL,    0, 0, "LEU ", 1, FALSE, 2},
                  {"Iso",      3, NULL,    0, 0, "ILE ", 1, FALSE, 2},
                  {"Ser",      3, NULL,    0, 0, "SER ", 1, FALSE, 2},
                  {"Thr",      3, NULL,    0, 0, "THR ", 1, FALSE, 2},
                  {"Cys",      3, NULL,    0, 0, "CYS ", 2, FALSE, 3},
                  {"Pro",      3, NULL,    0, 0, "PRO ", 1, FALSE, 2},
                  {"Phe",      3, NULL,    0, 0, "PHE ", 1, FALSE, 2},
                  {"Tyr",      3, NULL,    0, 0, "TYR ", 1, FALSE, 2},
                  {"Try",      3, NULL,    0, 0, "TRP ", 1, FALSE, 2},
                  {"His",      3, NULL,    0, 0, "HIS ", 2, FALSE, 3},
                  {"Aspartic", 8, NULL,    0, 0, "ASP ", 2, FALSE, 3},
                  {"Asparagi", 8, NULL,    0, 0, "ASN ", 1, FALSE, 2},
                  {"Glutamic", 8, NULL,    0, 0, "GLU ", 2, FALSE, 3},
                  {"Glutamin", 8, NULL,    0, 0, "GLN ", 1, FALSE, 2},
                  {"Methioni", 8, NULL,    0, 0, "MET ", 1, FALSE, 2},
                  {"Lys",      3, NULL,    0, 0, "LYS ", 1, FALSE, 2},
                  {"Arg",      3, NULL,    0, 0, "ARG ", 1, FALSE, 2},
                  {"Orn",      3, NULL,    0, 0, "ORN ", 1, TRUE,  2},
                  {"MethylAl", 8, NULL,    0, 0, "AIB ", 1, TRUE,  2},
                  {"Pyr",      3, NULL,    0, 0, "GLU ", 1, FALSE, 2},
                  {"Formyl",   6, NULL,    0, 0, "FOR ", 1, TRUE,  2},
                  {"Acetyl",   6, NULL,    0, 0, "ACE ", 1, TRUE,  2},
                  {"N-MeAmid", 8, NULL,    0, 0, "VAL ", 1, TRUE,  2},
                  {"N-Term",   6, "AIB",   3, 1, NULL,   2, TRUE,  3},
                  {"N-Term",   6, NULL,    0, 1, NULL,   2, FALSE, 3},
                  {"N-Term",   6, NULL,    0, 1, NULL,   3, FALSE, 4},
                  {"C-Term",   6, "Amide", 5, 0, "CTER", 2, FALSE, 3},
                  {"C-Term",   6, "AIB",   3, 1, NULL,   2, TRUE,  3},
                  {"C-Term",   6, "ORN",   3, 1, NULL,   2, TRUE,  3},
                  {"C-Term",   6, NULL,    0, 1, NULL,   2, FALSE, 3},
                  {"C-Term",   6, NULL,    0, 1, NULL,   3, FALSE, 4},
                  {"R-Aden",   6, NULL,    0, 0, "  A ", 1, FALSE, 2},
                  {"R-Guan",   6, NULL,    0, 0, "  G ", 1, FALSE, 2},
                  {"R-Cyto",   6, NULL,    0, 0, "  C ", 1, FALSE, 2},
                  {"R-Urac",   6, NULL,    0, 0, "  U ", 1, FALSE, 2},
                  {"D-Aden",   6, NULL,    0, 0, " DA ", 1, FALSE, 2},
                  {"D-Guan",   6, NULL,    0, 0, " DG ", 1, FALSE, 2},
                  {"D-Cyto",   6, NULL,    0, 0, " DC ", 1, FALSE, 2},
                  {"D-Urac",   6, NULL,    0, 0, " DU ", 1, FALSE, 2},
                  {"D-Thym",   6, NULL,    0, 0, " DT ", 1, FALSE, 2},
                  
                  {"R-Phos",   6, NULL,    0, 0, "PHO ", 1, TRUE,  2},
                  {"R-5'-Hyd", 8, NULL,    0, 0, "HYD ", 1, TRUE,  2},
                  {"R-5'-Pho", 8, NULL,    0, 0, "PHO ", 1, TRUE,  2},
                  {"R-3'-Hyd", 8, NULL,    0, 0, "HYD ", 1, TRUE,  2},
                  {"R-3'-Pho", 8, NULL,    0, 0, "PHO ", 1, TRUE,  2},
                  {"D-Phos",   6, NULL,    0, 0, "PHO ", 1, TRUE,  2},
                  {"D-5'-Hyd", 8, NULL,    0, 0, "HYD ", 1, TRUE,  2},
                  {"D-5'-Pho", 8, NULL,    0, 0, "PHO ", 1, TRUE,  2},
                  {"D-3'-Hyd", 8, NULL,    0, 0, "HYD ", 1, TRUE,  2},
                  {"D-3'-Pho", 8, NULL,    0, 0, "PHO ", 1, TRUE,  2},
                  
                  {"TIP3P",    5, NULL,    0, 0, "HOH ", 1, TRUE,  2},
                  {"Li+",      3, NULL,    0, 0, "LI  ", 0, TRUE,  3},
                  {"Na+",      3, NULL,    0, 0, "NA  ", 0, TRUE,  3},
                  {"K+",       2, NULL,    0, 0, "K   ", 0, TRUE,  3},
                  {"Rb+",      3, NULL,    0, 0, "RB  ", 0, TRUE,  3},
                  {"Cs+",      3, NULL,    0, 0, "CS  ", 0, TRUE,  3},
                  {"Mg+",      3, NULL,    0, 0, "MG  ", 0, TRUE,  3},
                  {"Ca+",      3, NULL,    0, 0, "CA  ", 0, TRUE,  3},
                  {"Zn+",      3, NULL,    0, 0, "ZN  ", 0, TRUE,  3},
                  {"Ba+",      3, NULL,    0, 0, "BA  ", 0, TRUE,  3},
                  {"Cl-",      3, NULL,    0, 0, "CL  ", 0, TRUE,  3},
                  {NULL,       0, NULL,    0, 0, NULL,   0, TRUE,  0}
   };

   int i;
   for(i=0; lookup[i].input1!=NULL; i++)
   {
      LOOKUP *l = &(lookup[i]);
      
      /* If we have the right number of fields                          */
      if(l->nfields == nwords)
      {
         /* If the first field matches                                  */
         if(!strncmp(l->input1, words[0], l->input1Len))
         {
            /* If we don't need to check the second field, or we do
               and it matches    
            */
            if((l->input2Len == 0) ||
               !strncmp(l->input2, words[1], l->input2Len))
            {
               /* We have a match.
                  If we have an output field specified, use that, 
                  otherwise use the one in the lookup structure
               */
               if(l->outField)
               {
                  strcpy(resnam, words[l->outField]);
               }
               else
               {
                  strcpy(resnam, l->output);
               }
               PADMINTERM(resnam, 4);
               

               /* Now get the atom field                                */
               if(!strncmp(words[l->atomField], "Oxygen", 6))
               {
                  strcpy(atnam, " O  ");
               }
               else if(!strncmp(words[l->atomField], "Hydrogen", 8))
               {
                  strcpy(atnam, " H  ");
               }
               else
               {
                  BOOL ion = FALSE;
                  char inputAtnam[MAXLABEL];
                  
                  strncpy(inputAtnam, words[l->atomField], MAXLABEL);

                  /* See if it's an ion                                 */
                  if(strchr(inputAtnam, '+') || strchr(inputAtnam, '-'))
                     ion = TRUE;
                  
                  /* Remove any charge information and up-case          */
                  TERMAT(inputAtnam, '+');
                  TERMAT(inputAtnam, '-');
                  UPPER(inputAtnam);

                  if(ion)    /* It's an ion                             */
                  {
                     /* If it's one character, we need a leading space  */
                     if(strlen(inputAtnam) == 1)
                     {
                        strcpy(atnam, " ");
                     }
                     strcat(atnam, inputAtnam);
                  }
                  else       /* It's a normal atom                      */
                  {
                     if(strlen(inputAtnam) == 4)
                     {
                        /* If it's 4 characters, move the last to the
                           start
                        */
                        atnam[0] = inputAtnam[3];
                        atnam[1] = '\0';
                        inputAtnam[3] = '\0';
                        strcat(atnam, inputAtnam);
                     }
                     else
                     {
                        /* Insert a leading space                       */
                        strcpy(atnam, " ");
                        strcat(atnam, inputAtnam);
                     }
                  }
                  
               }
               PADMINTERM(atnam, 4);

               return(TRUE);
            }
         }
      }
      
   }

   return(FALSE);
}




/************************************************************************/
void Usage(void)
{
   fprintf(stderr,"Usage: \n");
}



/************************************************************************/
PDB *ReadTinkerAsPDB(FILE *in, FILE *paramFp, char *header)
{
   char buffer[MAXBUFF];
   PDB  *pdb = NULL,
        *p   = NULL;
   char tinkerAtomTypesResnam[MAXATOMTYPES][MAXLABEL],
        tinkerAtomTypesAtnam[MAXATOMTYPES][MAXLABEL];
   BOOL tinkerAtomTypesIsHet[MAXATOMTYPES];
   int atnum, atomType;
   char tinkerAtnam[MAXLABEL];
   REAL x, y, z;
      

   ReadTinkerAtomTypes(paramFp, 
                       tinkerAtomTypesResnam, 
                       tinkerAtomTypesAtnam, 
                       tinkerAtomTypesIsHet, MAXATOMTYPES);

   
   /* Read the header line                                              */
   fgets(buffer, MAXBUFF, in);
   strcpy(header, buffer+8);

   /* Read the coordinates. We will put the atom type in formal_charge  */
   while(fgets(buffer, MAXBUFF, in))
   {
      if(pdb==NULL)
      {
         INIT(pdb, PDB);
         p=pdb;
      }
      else
      {
         ALLOCNEXT(p, PDB);
      }

      if(p==NULL)
      {
         FREELIST(pdb, PDB);
         return(FALSE);
      }
         
      fsscanf(buffer,"%6d%2x%2s%1x%12.6lf%12.6lf%12.6lf%1x%5d",
              &atnum, tinkerAtnam, &x, &y, &z, &atomType);

      PopulatePDBRecord(p, atnum, x, y, z, 
                        tinkerAtomTypesResnam[atomType],
                        tinkerAtomTypesAtnam[atomType],
                        tinkerAtomTypesIsHet[atomType]);
   }

   FixHydrogens(pdb);
   FixCterOxygens(pdb);
   FixAtomNames(pdb);
   FixILECD1(pdb);
   
   return(pdb);
}

void FixHydrogens(PDB *pdb)
{
   PDB *p, 
       *firstHydrogen;

   BOOL inHydrogens    = FALSE;
   int  hydrogenNumber = 2;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!inHydrogens)    /* Not currently in a block of hydrogens   */
      {
         if(p->atnam[0] == 'H')
         {
            inHydrogens    = TRUE;
            firstHydrogen  = p;
            hydrogenNumber = 2;
         }
      }
      else                /* Already in a block of hydrogens         */
      {
         if(p->atnam[0] != 'H')
         {
            /* Just come out of a block of hydrogens                 
               If this wasn't the atom immediately after the current
               firstHydrogen, then we need to put a '1' into the
               name of the firstHydrogen
            */
#ifdef DEBUG
            fprintf(stdout,"START\n");
            blWritePDBRecord(stdout,firstHydrogen);
            blWritePDBRecord(stdout,p);
            fprintf(stdout,"STOP\n");
#endif
            
            if(firstHydrogen->next != p)
            {
               InsertNumberInAtnam(firstHydrogen, 1);
            }
            
            inHydrogens    = FALSE;
            hydrogenNumber = 2;
         }
         else /* Still in hydrogens, but check if label has changed  */
         {
            if(strncmp(p->atnam, firstHydrogen->atnam, 4))
            {
               /* Label has changed
                  If this wasn't the atom immediately after the current
                  firstHydrogen, then we need to put a '1' into the
                  name of the firstHydrogen
               */
               if(firstHydrogen->next != p)
               {
                  InsertNumberInAtnam(firstHydrogen, 1);
               }
               /* Update the firstHydrogen to this atom since it's the
                  start of a new block
               */
               firstHydrogen  = p;
               hydrogenNumber = 2;
            }
            else  /* Label is the same                               */
            {
               /* We need to update the hydrogen atom label          */
               InsertNumberInAtnam(p, hydrogenNumber++);
            }
         }
      }
   }
}

/************************************************************************/
void InsertNumberInAtnam(PDB *p, int hydrogenNumber)
{
   char buffer[MAXLABEL],
        *ptr;

   sprintf(buffer, "%d", hydrogenNumber);

   if((ptr=strchr(p->atnam, ' '))!=NULL)
   {
      *ptr = buffer[0];
   }

   if((ptr=strchr(p->atnam_raw+1, ' '))!=NULL)
   {
      *ptr = buffer[0];
   }
   else if((ptr=strchr(p->atnam_raw, ' '))!=NULL)
   {
      *ptr = buffer[0];
   }
}



/************************************************************************/
void FixCterOxygens(PDB *pdb)
{
   PDB *p, *q, 
       *nextRes;

   for(p=pdb; p!=NULL; p=nextRes)
   {
      BOOL GotOXT = FALSE;
      
      nextRes = blFindNextResidue(p);
      for(q=p; q!=nextRes; NEXT(q))
      {
         if(!strncmp(q->atnam, "OXT ", 4))
         {
            if(!GotOXT)
            {
               strcpy(q->atnam,     "O   ");
               strcpy(q->atnam_raw, " O  ");
            }
            GotOXT = TRUE;
         }
      }
   }
}


/************************************************************************/
void FixILECD1(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->resnam, "ILE", 3) &&
         !strncmp(p->atnam,  "CD  ", 4))
      {
         strcpy(p->atnam,     "CD1 ");
         strcpy(p->atnam_raw, " CD1");
      }
   }
}


/************************************************************************/
void FixAtomNames(PDB *pdb)
{
   PDB *p, *q, 
       *nextRes;

   for(p=pdb; p!=NULL; p=nextRes)
   {
      nextRes = blFindNextResidue(p);

      if(!strncmp(p->resnam, "ASP", 3) ||
         !strncmp(p->resnam, "GLU", 3))
      {
         for(q=p; q!=nextRes; NEXT(q))
         {
            doFixAtomName(p, nextRes, "OD  ", 4);
            doFixAtomName(p, nextRes, "OE  ", 4);
         }
      }

      if(!strncmp(p->resnam, "TYR", 3) ||
         !strncmp(p->resnam, "PHE", 3))
      {
         doFixAtomName(p, nextRes, "CD  ", 4);
         doFixAtomName(p, nextRes, "CE  ", 4);
      }
      
      if(!strncmp(p->resnam, "ARG", 3))
      {
         doFixAtomName(p, nextRes, "NH  ", 4);
      }
      
   }
}

void doFixAtomName(PDB *start, PDB *stop, char *atom, int nChars)
{
   PDB *p;
   int count = 1;
   
   for(p=start; p!=stop; NEXT(p))
   {
      if(!strncmp(p->atnam, atom, nChars))
      {
         InsertNumberInAtnam(p, count++);
      }
   }
}



/************************************************************************/
void PopulatePDBRecord(PDB *p, int atnum, REAL x, REAL y, REAL z,
                       char *resnam, char *atnam, BOOL isHet)
{
   static int resnum = 0;
   CLEAR_PDB(p);
   strcpy(p->record_type, (isHet?"HETATM":"ATOM  "));
   p->atnum = atnum;
   strcpy(p->atnam_raw, atnam);
   strcpy(p->atnam, (atnam[0]==' '?atnam+1:atnam));
   PADMINTERM(p->atnam, 4);
   strcpy(p->resnam, resnam);
   PADMINTERM(p->resnam, 4);
   if(!strncmp(atnam, " N  ", 4))
      resnum++;
   p->resnum = resnum;
   p->x = x;   p->y = y;   p->z = z;
   p->occ = 1.0;
   blSetElementSymbolFromAtomName(p->element, atnam);
}


/************************************************************************/
/*>void DoChain(PDB *pdb, char **chains, BOOL BumpChainOnHet)
   ---------------------------------------------------------
*//**

   \param[in,out]  *pdb            PDB linked list
   \param[in]      *chains         Chain labels (or blank string)
   \param[in]      BumpChainOnHet  Bump the chain label when a HETATM
                                   is found

   Do the actual chain naming.

   *** CODE TAKEN FROM pdbchain.c ***

-  12.07.94 Original    By: ACRM
-  25.07.94 Only increments ch if *ch != \0
-  04.01.95 Added check on HETATM records
-  27.01.95 ChainNum count now mod 26 so labels will cycle A-Z
-  16.10.95 Handles BumpChainOnHet
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  05.03.15 Replaced blFindEndPDB() with blFindNextResidue()
-  10.03.15 Chains is now an array
*/
void DoChain(PDB *pdb, char **chains, BOOL BumpChainOnHet)
{
   PDB  *p,
        *start,
        *end,
        *LastStart = NULL,
        *N         = NULL,
        *C         = NULL,
        *CPrev     = NULL,
        *CAPrev    = NULL,
        *CA        = NULL;
   int  ChainNum   = 0,
        ChainIndex = 0;
   char chain[MAXCHAINLABEL];
   BOOL NewChain;
   

   if((chains!=NULL) && chains[ChainIndex][0])
      strcpy(chain, chains[ChainIndex++]);
   else
      strcpy(chain, "A");
   
   for(start=pdb; start!=NULL; start=end)
   {
      NewChain = FALSE;
      end = blFindNextResidue(start);

      CA = N = C = NULL;
      
      for(p=start; p!=end; NEXT(p))
      {
         if(!strncmp(p->atnam,"CA  ",4)) CA  = p;
         if(!strncmp(p->atnam,"N   ",4)) N   = p;
         if(!strncmp(p->atnam,"C   ",4)) C   = p;
      }

      if(CPrev != NULL && N != NULL)
      {
         /* A C was defined in the last residue and an N in this one
            Calc C-N distance
         */
         if(DISTSQ(CPrev, N) > CNDISTSQ)
            NewChain = TRUE;
      }
      else if(CAPrev != NULL && CA != NULL)
      {
         /* No C-N connection, but a CAs found
            Calc CA-CA distance
         */
         if(DISTSQ(CAPrev, CA) > CADISTSQ)
            NewChain = TRUE;
      }
      else if(LastStart != NULL)
      {
         char buffer[80],
              atoms[80];
         
         /* Build string specifying faulty residues                     */
         if((CPrev == NULL || CAPrev == NULL) &&
            (N     == NULL || CA     == NULL))
            sprintf(buffer,"residues %s.%d%c and %s.%d%c",
                    LastStart->chain,LastStart->resnum,
                    LastStart->insert[0],
                    start->chain,start->resnum,start->insert[0]);
         else if(CPrev == NULL || CAPrev == NULL)
            sprintf(buffer,"residue %s.%d%c",
                    LastStart->chain,LastStart->resnum,
                    LastStart->insert[0]);
         else
            sprintf(buffer,"residue %s.%d%c",
                    start->chain,start->resnum,start->insert[0]);

         /* Build string specifying faulty atoms                        */
         atoms[0] = '\0';
         if(CAPrev == NULL || CA == NULL) strcat(atoms,"CA ");
         if(N      == NULL)               strcat(atoms,"N ");
         if(CPrev  == NULL)               strcat(atoms,"C ");

         /* Print warning message                                       */
         if((LastStart != NULL && 
             strncmp(LastStart->record_type, "HETATM", 6)) &&
            (start     != NULL && 
             strncmp(start->record_type,     "HETATM", 6)))
            fprintf(stderr, "Warning: Atoms missing in %s: %s\n",
                    buffer, atoms);

         if(BumpChainOnHet &&
            (LastStart != NULL && 
             !strncmp(LastStart->record_type, "HETATM", 6)) &&
            (start     != NULL && 
             !strncmp(start->record_type,     "ATOM  ",6)))
            NewChain = TRUE;
      }

      /* If we've changed chain, set the new chain name                 */
      if(NewChain)
      {
         ChainNum++;
         
         if((chains!=NULL) && chains[ChainIndex][0])
         {
            strcpy(chain,chains[ChainIndex++]);
         }
         else
         {
            strcpy(chain, GetChainLabel(ChainNum));
         }
      }

      /* Copy the name into this residue                                */
      for(p=start; p!=end; NEXT(p))
         strcpy(p->chain, chain);
      
      /* Set pointers for next residue                                  */
      CAPrev    = CA;
      CPrev     = C;
      LastStart = start;
   }
}


/************************************************************************/
/*>char *GetChainLabel(int ChainNum)
   ---------------------------------
*//**
   \param[in]  ChainNum    Chain number
   \return                 Chain label 

   Converts a chain number (>=0) into a chain label. Chain labels run
   from A-Z, a-z, 1-9, 0, and then 63 onwards as multi-character strings

   *** CODE TAKEN FROM pdbchain.c ***

-  10.03.15 Original   By: ACRM
*/
char *GetChainLabel(int ChainNum)
{
   static char chain[MAXCHAINLABEL];
   
   if(ChainNum < 26)
   {
      chain[0] = (char)(65 + ChainNum);
      chain[1] = '\0';
   }
   else if(ChainNum < 52)
   {
      chain[0] = (char)(97 + (ChainNum-26));
      chain[1] = '\0';
   }
   else if(ChainNum < 61)
   {
      sprintf(chain,"%d", ChainNum-51);
   }
   else if(ChainNum == 61)
   {
      strcpy(chain,"0");
   }
   else
   {
      sprintf(chain,"%d", ChainNum);
   }
   
   return(chain);
}

/************************************************************************/
void RenumberResidues(PDB *pdb)
{
   PDB  *p;
   int  resnum   = 0,
        LastRes;
   char LastInsert[MAXLABEL],
        LastChain[MAXCHAINLABEL];

   LastRes       = (-9999);
   LastInsert[0] = '\0';
   LastChain[0]  = '\0';
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* Increment resnum if we have changed residue                    */
      if((p->resnum != LastRes) ||
         !INSERTMATCH(p->insert, LastInsert))
      {
         LastRes = p->resnum;
         strcpy(LastInsert, p->insert);
         
         resnum++;
      }

      /* See if we've changed chain                                     */
      if(!CHAINMATCH(p->chain, LastChain))
      {
         resnum = 1;
         strcpy(LastChain, p->chain);
      }
      
      /* Set the residue number                                         */
      p->resnum = resnum;

      /* Set the insert code to a blank                                 */
      strcpy(p->insert, " ");
   }
}
