#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include "fasta.h"

int processFasta(char *filename, double *timeTaken, int *arraySize)
{
	FILE *fp;
	FASTArecord fRecord;
	int recordNumber = 0, status;
	int eofSeen = 0;
	clock_t startTime, endTime;
	FASTArecord *recordArray = NULL;
	// int arraySize = 1000; // Initial array size

	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "Failure opening %s : %s\n", filename, strerror(errno));
		return -1;
	}

	/** record the time now, before we do the work */
	startTime = clock();
	int count = 0;

	/*allocating memory for an array of FASTA records*/
	recordArray = (FASTArecord *)malloc(*arraySize * sizeof(FASTArecord));

	do
	{
		/** print a '.' every 10,000 records so
		 * we know something is happening */
		if ((recordNumber % 10000) == 0)
		{
			printf(".");
			fflush(stdout);
		}

		fastaInitializeRecord(&fRecord);

		status = fastaReadRecord(fp, &fRecord);
		if (status == 0)
		{
			eofSeen = 1;
		}
		else if (status > 0)
		{
			if (recordNumber == *arraySize)
			{
				// doubling array size and reallocating memory
				*arraySize *= 2;
				recordArray = (FASTArecord *)realloc(recordArray, *arraySize * sizeof(FASTArecord));
			}
			// Copying into arrray
			recordArray[recordNumber] = fRecord;
			recordNumber++;
			count++;
		}
		else
		{
			fprintf(stderr, "Error: failure at line %d of '%s'\n", recordNumber, filename);
			return -1;
		}

	} while (!eofSeen);

	endTime = clock();
	(*timeTaken) = ((double)(endTime - startTime)) / CLOCKS_PER_SEC;

	fclose(fp);

	// freeing memory for descriptions and sequences in recordArray
	for (int i = 0; i < recordNumber; i++)
	{
		free(recordArray[i].description);
		free(recordArray[i].sequence);
	}

	// freeing recordArray
	free(recordArray);
	printf("........... %d FASTA records", count);
	return recordNumber;
}

int processFastaRepeatedly(char *filename, long repeatsRequested)
{
	// double timeThisIterationInSeconds;
	double totalTimeInSeconds = 0;
	int minutesPortion;
	int status;
	long i;
	int arraySize = 1000;

	for (i = 0; i < repeatsRequested; i++)
	{
		double timeThisIterationInSeconds;
		status = processFasta(filename, &timeThisIterationInSeconds, &arraySize);
		if (status < 0)
		{
			return -1;
		}
		// calculating memory waste
		int memoryWaste = ((arraySize - status) * 100) / arraySize;

		printf(" -- %d allocated (%d%% waste)\n", arraySize, memoryWaste);

		totalTimeInSeconds += timeThisIterationInSeconds;
	}

	printf("%.3f seconds taken for processing total\n", totalTimeInSeconds);

	totalTimeInSeconds /= (double)repeatsRequested;

	minutesPortion = (int)(totalTimeInSeconds / 60);
	totalTimeInSeconds = totalTimeInSeconds - (60 * minutesPortion);
	printf("On average: %d minutes, %.3f second per run\n", minutesPortion, totalTimeInSeconds);

	return status;
}

void usage(char *progname)
{
	fprintf(stderr, "%s [<OPTIONS>] <file> [ <file> ...]\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Prints timing of loading and storing FASTA records.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Options: \n");
	fprintf(stderr, "-R <REPEATS> : Number of times to repeat load.\n");
	fprintf(stderr, "             : Time reported will be the average time.\n");
	fprintf(stderr, "\n");
}

int main(int argc, char **argv)
{
	int i, recordsProcessed = 0;
	long repeatsRequested = 1;

	for (i = 1; i < argc; i++)
	{
		if (argv[i][0] == '-')
		{
			if (argv[i][1] == 'R')
			{
				if (i >= argc)
				{
					fprintf(stderr, "Error: need an argument for repeats requested\n");
					return 1;
				}
				if (sscanf(argv[++i], "%ld", &repeatsRequested) != 1)
				{
					fprintf(stderr, "Error: cannot parse repeats requested from '%s'\n", argv[i]);
					return 1;
				}
			}
			else
			{
				fprintf(stderr, "Error: unknown option '%s'\n", argv[i]);
				usage(argv[0]);
			}
		}
		else
		{
			recordsProcessed = processFastaRepeatedly(argv[i], repeatsRequested);
			if (recordsProcessed < 0)
			{
				fprintf(stderr, "Error: Processing '%s' failed -- exiting\n", argv[i]);
				return 1;
			}
			/*printf("%d records processed from '%s'\n",
					recordsProcessed, argv[i]);*/
		}
	}

	if (recordsProcessed == 0)
	{
		fprintf(stderr, "No data processed -- provide the name of a file on the command line\n");
		usage(argv[0]);
		return 1;
	}

	return 0;
}