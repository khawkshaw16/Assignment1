#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <math.h>
using namespace std;
//set up the bigram calc for the nucleotides
void bigramCalc(char nuc, char previous, int& AA, int& AT, int& AC, int& AG, int& TT, 
	int& TA, int& TC, int& TG, int& CC, int& CA, int& CT, int& CG, int& GG, int& GA, 
	int& GC, int& GT) {
	if (previous == '\0') {
		return;
	}
	else if (nuc == 'A' && previous == 'A') {
		AA++;
	}
	else if (nuc == 'A' && previous == 'T') {
		AT++;
	}
	else if (nuc == 'A' && previous == 'C') {
		AC++;
	}
	else if (nuc == 'A' && previous == 'G') {
		AG++;
	}
	else if (nuc == 'T' && previous == 'T') {
		TT++;
	}
	else if (nuc == 'T' && previous == 'A') {
		TA++;
	}
	else if (nuc == 'T' && previous == 'C') {
		TC++;
	}
	else if (nuc == 'T' && previous == 'G') {
		TG++;
	}
	else if (nuc == 'C' && previous == 'C') {
		CC++;
	}
	else if (nuc == 'C' && previous == 'A') {
		CA++;
	}
	else if (nuc == 'C' && previous == 'T') {
		CT++;
	}
	else if (nuc == 'C' && previous == 'G') {
		CG++;
	}
	else if (nuc == 'G' && previous == 'G') {
		GG++;
	}
	else if (nuc == 'G' && previous == 'A') {
		GA++;
	}
	else if (nuc == 'G' && previous == 'C') {
		GC++;
	}
	else if (nuc == 'G' && previous == 'T') {
		GT++;
	}

}

//gaussian distribution

double gaussianDis(double& variance, double& currentMean) {
	double a = ((double) rand() / (RAND_MAX));
	double b = ((double)rand() / (RAND_MAX));

	double c = sqrt(-2 * (log(a) / log(exp(1)) * cos(2 * M_PI * b)));

	double d = variance * c + currentMean;
	return d;
}
//main function
int main() {
	//set up methods to find mean, sum, variance, stdv
	int lengthTotal=0, aTotal=0, cTotal=0, tTotal=0, gTotal=0, AA=0, AT=0, AC=0, AG=0, TT=0, TA=0, TC=0, TG=0, CC=0, CA=0, CT=0, CG=0, GG=0, GA=0, GC=0, GT=0;
	double meanDivisor = 0;
	int currentLineLength = 0;
	double currentMean = 0;
	double previousMean = 0;
	int numLines = 0;
	double variance = 0;
	
	double stdv = 0;
	

	//input file
	ifstream infile("dna.txt");
	//output file
	ofstream outfile("kellyhawkshaw.out");
	char nuc;
	char previous = '\0';
	//while loop to set up string list with the sequence of nucleotides
		while (infile) {
			char nuc;
			
				infile >> nuc;
				
				nuc = toupper(nuc);

				if (nuc == 'A') {
					currentLineLength++;
					aTotal++;
					lengthTotal++;
					bigramCalc(nuc, previous, AA, AT, AC, AG, TT, TA, TC, TG, CC, CA, CT, CG, GG, GA, GC, GT);
				}
				else if (nuc == 'C') {
					currentLineLength++;
					cTotal++;
					lengthTotal++;
					bigramCalc(nuc, previous, AA, AT, AC, AG, TT, TA, TC, TG, CC, CA, CT, CG, GG, GA, GC, GT);
				}

				else if (nuc == 'T') {
					currentLineLength++;
					tTotal++;
					lengthTotal++;
					bigramCalc(nuc, previous, AA, AT, AC, AG, TT, TA, TC, TG, CC, CA, CT, CG, GG, GA, GC, GT);
				}

				else if (nuc == 'G') {
					currentLineLength++;
					gTotal++;
					lengthTotal++;
					bigramCalc(nuc, previous, AA, AT, AC, AG, TT, TA, TC, TG, CC, CA, CT, CG, GG, GA, GC, GT);
				}
				else if (nuc == '\n') {
					//mean
					meanDivisor++;
					currentMean = (previousMean + ((currentLineLength - previousMean) / meanDivisor));


				}

				currentLineLength = 0;
				previous = nuc;


			}
			
		//variance 
		variance = lengthTotal / meanDivisor - (currentMean * currentMean);
		//standard deviation
		stdv = sqrt(variance);


		//output for the values
		outfile << "Kelly Hawkshaw, ID:2328274" << endl;
		outfile << "Mean Line Length: " << currentMean << endl;
		outfile << "Total Number of Nucleotides: " << lengthTotal << endl;
		outfile << "Standard Deviation : " << stdv << endl;
		outfile << "Variance: " << variance << endl;
		outfile << "Bigram AA Frequency : " << ((double) AA)/lengthTotal << endl;
		outfile << "Bigram AT Frequency : " << ((double)AT) / lengthTotal << endl;
		outfile << "Bigram AC Frequency : " << ((double)AC) / lengthTotal << endl;
		outfile << "Bigram AG Frequency : " << ((double)AG) / lengthTotal << endl;
		outfile << "Bigram TT Frequency : " << ((double)TT) / lengthTotal << endl;
		outfile << "Bigram TA Frequency : " << ((double)TA) / lengthTotal << endl;
		outfile << "Bigram TC Frequency : " << ((double)TG) / lengthTotal << endl;
		outfile << "Bigram CC Frequency : " << ((double)CC) / lengthTotal << endl;
		outfile << "Bigram CA Frequency : " << ((double)CA) / lengthTotal << endl;
		outfile << "Bigram CT Frequency : " << ((double)CT) / lengthTotal << endl;
		outfile << "Bigram CG Frequency : " << ((double)CG) / lengthTotal << endl;
		outfile << "Bigram GG Frequency : " << ((double)GG) / lengthTotal << endl;
		outfile << "Bigram GA Frequency : " << ((double)GA) / lengthTotal << endl;
		outfile << "Bigram GC Frequency : " << ((double)GC) / lengthTotal << endl;
		outfile << "Bigram GT Frequency : " << ((double)GT) / lengthTotal << endl;
		
		//Generating 1000 DNA strings using Gaussian Distribution
		int userTok;


			for (int i = 0; i < 1000; ++i) {
				double x = gaussianDis(variance, currentMean);
				for (int j = 0; j < x; ++j) {
					int y = (rand() % 4) + 1;
					if (y == 1) {
						outfile << 'A';
					}
					else if (y == 2) {
						outfile << 'T';
					}
					else if (y == 3) {
						outfile << 'C';
					}
					else if (y == 4) {
						outfile << 'G';
					}
				}
				outfile << '\n';
			}
			//Ask if user wants to process another list
			cout << "Would you like to process another list? Any int for yes, -1 for no" << endl;
			cin >> userTok;

			while (userTok != -1) {
				for (int i = 0; i < 1000; ++i) {
					double x = gaussianDis(variance, currentMean);
					for (int j = 0; j < x; ++j) {
						int y = (rand() % 4) + 1;
						if (y == 1) {
							outfile << 'A';
						}
						else if (y == 2) {
							outfile << 'T';
						}
						else if (y == 3) {
							outfile << 'C';
						}
						else if (y == 4) {
							outfile << 'G';
						}
					}
					outfile << '\n';
				}
				//Ask if user wants to process another list
				cout << "Would you like to process another list? Any int for yes, -1 for no" << endl;
				cin >> userTok;
			}


	//closing input file
	infile.close();

	return 0;
}

