/*************************************************************************

   Program:    fixoverlap
   File:       fixoverlap.c
   
   Version:    V1.0
   Date:       19.12.19
   Function:   Checks a Tinker xyz file for hydrogen atoms with identical
               coordinates. If found moves the second one slightly
   
   Copyright:  (c) UCL / Prof. Andrew C. R. Martin 2019
   Author:     Prof. Andrew C. R. Martin
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
   V1.0   19.12.19  Original   By: ACRM

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
#define MAXWORD         16
#define SMALL      0.00001
typedef struct _tinkerxyz
{
   struct _tinkerxyz *next;
   REAL x,y,z;
   int  atnum,
        type,
        connect[8];
   char atnam[MAXLABEL];
}  TINKERXYZ;


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
TINKERXYZ *ReadTinkerXYZ(FILE *fp, int *natoms);
void WriteTinkerXYZ(FILE *fp, int natoms, TINKERXYZ *xyz);
void FixOverlaps(TINKERXYZ *xyz);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   FILE *in      = stdin,
        *out     = stdout;
   int  natoms;
   TINKERXYZ *xyz = NULL;
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         
         if((xyz=ReadTinkerXYZ(in, &natoms))==NULL)
         {
            fprintf(stderr,"Error: No atoms read from Tinker XYZ \
file\n");
            return(1);
         }

         FixOverlaps(xyz);
         
         WriteTinkerXYZ(out, natoms, xyz);
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
void FixOverlaps(TINKERXYZ *xyz)
{
   TINKERXYZ *a, *b;
   for(a=xyz; a!=NULL; NEXT(a))
   {
      for(b=a->next; b!=NULL; NEXT(b))
      {
         if((ABS(a->x - b->x) < SMALL) &&
            (ABS(a->y - b->y) < SMALL) &&
            (ABS(a->z - b->z) < SMALL))
         {
            fprintf(stderr, "Fixing %d\n", b->atnum);
            b->x += 1.0;
         }
      }
   }
}


/************************************************************************/
void WriteTinkerXYZ(FILE *fp, int natoms, TINKERXYZ *xyz)
{
   TINKERXYZ *t;
   int i;
   
   
   fprintf(fp, "%6d\n", natoms);
   for(t=xyz; t!=NULL; NEXT(t))
   {
      fprintf(fp, "%6d  %-3s%12.6f%12.6f%12.6f%6d",
              t->atnum, t->atnam,
              t->x, t->y, t->z,
              t->type);
      for(i=0; i<4; i++)
      {
         if(!t->connect[i])
            break;
         
         fprintf(fp, "%6d", t->connect[i]);
      }
      fprintf(fp, "\n");
   }
}


/************************************************************************/
TINKERXYZ *ReadTinkerXYZ(FILE *fp, int *natoms)
{
   TINKERXYZ *xyz = NULL,
             *t   = NULL;
   char      buffer[MAXBUFF],
      *chp,
      word[MAXWORD];
   
   int i;
   
   fgets(buffer, MAXBUFF, fp);
   if(!sscanf(buffer, "%d", natoms))
      return(NULL);

   if(natoms)
   {
      while(fgets(buffer, MAXBUFF, fp))
      {
         TERMINATE(buffer);
         
         if(xyz==NULL)
         {
            INIT(xyz, TINKERXYZ);
            t = xyz;
         }
         else
         {
            ALLOCNEXT(t, TINKERXYZ);
         }
         if(t==NULL)
         {
            FREELIST(xyz, TINKERXYZ);
            return(NULL);
         }

         for(i=0; i<4; i++)
         {
            t->connect[i] = 0;
         }

         sscanf(buffer, "%d%s%lf%lf%lf%d",
                &(t->atnum), t->atnam,
                &(t->x), &(t->y), &(t->z),
                &(t->type));

         chp = buffer+53; /* Start of connects */
         
         i=0;
         do
         {
            chp = blGetWord(chp, word, MAXWORD);
            sscanf(word, "%d", &(t->connect[i++]));
         } while(chp != NULL);
         
      }
   }

   return(xyz);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*/
void Usage(void)
{
   
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, 
                     char *infile, char *outfile)
   --------------------------------------------------------------
   Input:   int    argc          Argument count
            char   **argv        Argument array
   Output:  char   *infile       Input file (or blank string)
            char   *outfile      Output file (or blank string)
   Returns: BOOL                 Success?

   Parse the command line

   19.12.19  Original   By: ACRM  
*/
BOOL ParseCmdLine(int argc, char **argv, 
                  char *infile, char *outfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
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
         /* Check that there are 0-2 arguments left                     */
         if(argc > 2)
            return(FALSE);

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


