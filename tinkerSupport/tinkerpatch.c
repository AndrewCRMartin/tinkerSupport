/*************************************************************************

   Program:    tinkerpatch
   File:       tinkerpatch
   
   Version:    V1.0
   Date:       17.09.15
   Function:   Patch the numbering in a Tinker PDB file with the numbering
               from the original PDB file
   
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
#define MAXLABEL         8

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *paramFile, 
                  char *infile, char *outfile);
void Usage(void);
BOOL tinkerpatch(PDB *pdbNew, PDB *pdbOld);
void FixResidueNames(PDB *pdb);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*/
int main(int argc, char **argv)
{
   char origFile[MAXBUFF],
        infile[MAXBUFF],
        outfile[MAXBUFF];
   FILE *in      = stdin,
        *out     = stdout,
        *fp      = NULL;
   PDB  *pdbOrig = NULL,
        *pdbNew  = NULL;
   int  natoms;
    
   if(ParseCmdLine(argc, argv, origFile, infile, outfile))
   {
      if((fp=fopen(origFile, "r"))==NULL)
      {
         fprintf(stderr,"Error: Unable to open original PDB \
file: %s\n", origFile);
         
         return(1);
      }
      
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdbOrig=blReadPDB(fp, &natoms))==NULL)
         {
            fprintf(stderr,"Error: No atoms read from original PDB \
file\n");
            return(1);
         }

         if((pdbNew=blReadPDB(in, &natoms))==NULL)
         {
            fprintf(stderr,"Error: No atoms read from minimized PDB \
file\n");
            return(1);
         }
         
            
         if(!tinkerpatch(pdbNew, pdbOrig))
         {
            fprintf(stderr,"Error: Patching failed\n");
            return(1);
         }

         blWritePDB(out, pdbNew);
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
/*>BOOL tinkerpatch(PDB *pdbNew, PDB *pdbOld)
   ------------------------------------------
*/
BOOL tinkerpatch(PDB *pdbNew, PDB *pdbOld)
{
   PDB *p,
       *pNewStart,
       *pNewStop,
       *pOldStart;
   int atnum = 0;
   
   FixResidueNames(pdbNew);
   
   pOldStart = pdbOld;
   for(pNewStart=pdbNew; pNewStart!=NULL; pNewStart=pNewStop)
   {
      if(pOldStart == NULL)
      {
         fprintf(stderr,"Error: Original structure ran out of \
residues!\n");
      }
      
      pNewStop = blFindNextResidue(pNewStart);
      for(p=pNewStart; p!=pNewStop; NEXT(p))
      {
         atnum++;
         
         if(strncmp(p->resnam, pOldStart->resnam, 4))
         {
            fprintf(stderr,"Error: residue names don't match!\n");
            blWritePDBRecord(stderr, p);
            blWritePDBRecord(stderr, pOldStart);
            return(FALSE);
         }
            
         strcpy(p->chain, pOldStart->chain);
         p->resnum = pOldStart->resnum;
         strcpy(p->insert, pOldStart->insert);
      }
      pOldStart = blFindNextResidue(pOldStart);
   }

   return(TRUE);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*/
void Usage(void)
{
   
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *origFile, 
                     char *infile, char *outfile)
   --------------------------------------------------------------
   Input:   int    argc          Argument count
            char   **argv        Argument array
   Output:  char   *origFile
            char   *infile       Input file (or blank string)
            char   *outfile      Output file (or blank string)
   Returns: BOOL                 Success?

   Parse the command line

   17.09.15  Original   By: ACRM  
*/
BOOL ParseCmdLine(int argc, char **argv, char *origFile, 
                  char *infile, char *outfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = origFile[0] = '\0';
   
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
         
         /* Copy the first to origFile                                 */
         strcpy(origFile, argv[0]);
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
/*>void FixResidueNames(PDB *pdb)
   ------------------------------
*/
void FixResidueNames(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->resnam, "CYX", 3))
         strcpy(p->resnam, "CYS ");
      p->occ = 1.0;
   }
}


