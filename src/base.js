/************************************************************
* Filename: base.js
* Abstract: Definition of base class for KaKs methods.

* Version: js-2.0
* Author: Fanix
* Date: December.11, 2022

* Adapted from: KaKs_Calculator/base.cpp and KaKs_Calculator/base.h
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.

* Modified Version: 2.0.2
* Modified Author: Kristian K Ullrich
* Modified Date: October.06, 2022

* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

*************************************************************/

VERSION = "2.0"
DNASIZE = 4
XSIZE = DNASIZE * DNASIZE
CODONLENGTH = 3
DNASIZE = 4
XSIZE = DNASIZE * DNASIZE
CODON = 64
NULL = 0
NA = -1
NCODE = 23
NNCODE = NCODE * 2
SMALLVALUE = 1e-6
NUMBER_OF_RATES = 6
MODELCOUNT = 14

function gammap(x, alpha) { return alpha * (1 - pow(x, -1 / alpha)) }
function square(x) { return x * x }
function min2(a, b) { return a < b ? a : b }
function max2(a, b) { return a > b ? a : b }
function SIGN(a, b) { return (b) >= 0.0 ? Math.abs(a) : -Math.abs(a) }


/******** Global variables ********/


/*						The Genetic Codes 
http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
	Last update of the Genetic Codes: October 06, 2022 */
let genetic_code = 1; //from 1 to 33
/* Genetic standard codon table, !=stop codon */

let transl_table = [
	"FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "1-Standard Code",
	"FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS!!VVVVAAAADDEEGGGG", "2-Vertebrate Mitochondrial Code",
	"FFLLSSSSYY!!CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "3-Yeast Mitochondrial Code",
	"FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "4-Mold Mitochondrial Code",
	"FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "5-Invertebrate Mitochondrial Code",
	"FFLLSSSSYYQQCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "6-Ciliate, Dasycladacean and Hexamita Code",
	"", "7-",
	"", "8-",
	"FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "9-Echinoderm and Flatworm Mitochondrial Code",
	"FFLLSSSSYY!!CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "10-Euplotid Nuclear Code",
	"FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "11-Bacterial and Plant Plastid Code",
	"FFLLSSSSYY!!CC!WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "12-Alternative Yeast Nuclear Code",
	"FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "13-Ascidian Mitochondrial Code",
	"FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "14-Alternative Flatworm Mitochondrial Code",
	"FFLLSSSSYY!QCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "15-Blepharisma Nuclear Code",
	"FFLLSSSSYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "16-Chlorophycean Mitochondrial Code",
	"", "17-",
	"", "18-",
	"", "19-",
	"", "20-",
	"FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "21-Trematode Mitochondrial Code",
	"FFLLSS!SYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "22-Scenedesmus obliquus mitochondrial Code",
	"FF!LSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "23-Thraustochytrium Mitochondrial Code",
	"FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "24-Rhabdopleuridae Mitochondrial Code",
	"FFLLSSSSYY!!CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "25-Candidate Division SR1 and Gracilibacteria Code",
	"FFLLSSSSYY!!CC!WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "26-Pachysolen tannophilus Nuclear Code",
	"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "27-Karyorelict Nuclear Code",
	"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "28-Condylostoma Nuclear Code",
	"FFLLSSSSYYYYCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "29-Mesodinium Nuclear Code",
	"FFLLSSSSYYEECC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "30-Peritrich Nuclear Code",
	"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "31-Blastocrithidia Nuclear Code",
	"", "32-",
	"FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "33-Cephalodiscidae Mitochondrial UAA-Tyr Code"
]
let seq_name			//sequences' name
let length = 0			//sequences' length
let GC = new Array(4)	//GC Content
//********End of Global variables**********

//Constructor function
class Base {
	constructor() {
		this.Si = new Array(5).fill(0);
		this.Vi = new Array(5).fill(0);
		this.L = new Array(5).fill(0);
		this.K = new Array(5);
		this.kappa = this.kappatc = this.kappaag = 0;
		this.MLResult = {
			result: "",
			AICc: 0,
			freq: new Array(CODON),
			rate: new Array(NUMBER_OF_RATES),
			w: 0,
			t: 0
		}

		this.t = 0;

		this.KAPPA = new Array(NUMBER_OF_RATES).fill(1)

		this.SEKa = this.SEKs = this.AICc = this.lnL = this.AkaikeWeight = NA;
		this.Ka = this.Ks = this.Sd = this.Nd = this.S = this.N = this.snp = this.t = this.kappa = 0;
		this.model = ""
	}

	/**********************************************************************
	* Function: getAminoAcid
	* Input Parameter: codon or codon's id
	* Output: Calculate the amino acid according to codon or codon's id.
	* Return Value: char 
	***********************************************************************/
	getAminoAcid(codon) {
		if (codon == +codon)
			return transl_table[2 * (genetic_code - 1)][codon]
		return transl_table[2 * (genetic_code - 1)][this.getID(codon)];
	}

	/**********************************
	* Function: getNumNonsense
	* Input Parameter: int
	* Output: get the number of nonsense codons
	* Return Value: int
	***********************************/
	getNumNonsense() {
		let num, i;
		for (num = i = 0; i < CODON; i++) {
			if (this.getAminoAcid(i) == '!') num++;
		}
		return num;
	}
	/********************************************
	* Function: getID
	* Input Parameter: codon
	* Output: Get codon's id in array of codon_table.
	* Return Value: int
	*********************************************/
	getID(codon) {
		return (this.convertChar(codon[0]) * XSIZE + this.convertChar(codon[1]) * DNASIZE + this.convertChar(codon[2]));
	}
	/********************************************
	* Function: getCodon
	* Input Parameter: int
	* Output: Get the codon according to id;
			  a reverse funtion of getID.
	* Return Value: string
	*********************************************/
	getCodon(IDcodon) {
		let codon = "TTT";
		if (IDcodon >= 0 && IDcodon < 64) {
			codon[0] = this.convertInt(IDcodon / 16);
			codon[1] = this.convertInt((IDcodon % 16) / 4);
			codon[2] = this.convertInt(IDcodon % 4);
		}
		return codon;
	}
	/*********************************************
	* Function: convertChar
	* Input Parameter: ch as char
	* Output: Convert a char-T,C,A,G into a digit
	*		  0,1,2,3, respectively.
	* Return Value: int.
	**********************************************/
	convertChar(ch) {
		let ret = -1;
		switch (ch) {
			case 'T': case 'U':
				ret = 0;
				break;
			case 'C':
				ret = 1;
				break;
			case 'A':
				ret = 2;
				break;
			case 'G':
				ret = 3;
				break;
		}
		return ret;
	}

	/********************************************
	* Function: convertInt
	* Input Parameter: int
	* Output: Convert a digit- 0,1,2,3 into a 
	*		  char-T,C,A,G, respectively.
	* Return Value: char
	*********************************************/
	convertInt(i) {
		let ch = '-';
		switch (i) {
			case 0:
				ch = 'T';
				break;
			case 1:
				ch = 'C';
				break;
			case 2:
				ch = 'A';
				break;
			case 3:
				ch = 'G';
				break;
		}
		return ch;
	}
	/********************************************
	* Function: stringtoUpper
	* Input Parameter: string
	* Output: upper string
	* Return Value: string
	*********************************************/

	stringtoUpper(str) {
		return str.toUpperCase();
	}
	/********************************************
	* Function: getRandom
	* Input Parameter: void
	* Output: Generate a radnom integer
	* Return Value: int
	*********************************************/
	getRandom() {
		return Math.random();
	}

	/********************************************
	* Function: initArray
	* Input Parameter: array of int/double, int, int/double(default=0)
	* Output: Init the array x[0...n-1]=value
	* Return Value: int
	*********************************************/
	initArray(x, n, value = 0) {
		let i;
		for (i = 0; i < n; i++) x[i] = value;
		return 0;
	}
	/********************************************
	* Function: sumArray
	* Input Parameter: double/int, int, int(default=0)
	* Output: Sum of array x[]
	* Return Value: double/int
	*********************************************/
	sumArray(x, end, begin = 0) {
		let i, sum = 0.;
		for (i = begin; i < end; sum += x[i], i++);
		return sum;
	}

	/********************************************
	* Function: norm
	* Input Parameter: array of double, int
	* Output: Sqrt of the sum of the elements' square 
			   sqrt(x0*x0 + x1*x1 + ...)
	* Return Value: double
	*********************************************/
	norm(x, n) {
		let i, t = 0;
		for (i = 0; i < n; t += square(x[i]), i++);
		return Math.sqrt(t);
	}

	/********************************************
	* Function: scaleArray
	* Input Parameter: double, array of double, int
	* Output: Elements in array are mutipled by scale 
	* Return Value: int
	*********************************************/
	scaleArray(scale, x, n) {
		let i;
		for (i = 0; i < n; i++) x[i] *= scale;
		return 1;
	}

	/********************************************
	* Function: innerp
	* Input Parameter: array, array, int
	* Output: Sum of 'n' products multiplied by 
				two elements in x[], y[].
	* Return Value: int
	*********************************************/
	innerp(x, y, n) {
		let i, t = 0;
		for (i = 0; i < n; t += x[i] * y[i], i++);
		return t;
	}

	/********************************************
	* Function: initIdentityMatrix
	* Input Parameter: array of double, int
	* Output: Set x[i,j]=0 when x!=j and 
				  x[i,j]=1 when x=j 
	* Return Value: int
	*********************************************/
	initIdentityMatrix(x, n) {
		let i, j
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; x[i * n + j] = 0, j++);
			x[i * n + i] = 1;
		}
		return x;
	}

	/*****************************************************
	* Function: parseOutput
	* Input Parameter: void
	* Output: Parse estimated results for outputing
	* Return Value: string
	
	  Order: "Sequence", "Method", "Ka", "Ks", "Ka/Ks", 
			 "P-Value(Fisher)", "Length", "S-Sites", "N-Sites", "Fold-Sites(0:2:4)",
			 "Substitutions", "S-Substitutions", "N-Substitutions", "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)", 
			 "Divergence-Time", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
			 "Model"
	******************************************************/
	parseOutput() {
		let result = [];
		//Sequence name
		result.push(seq_name)
		//Method name
		result.push(this.name)
		//Ka
		result.push(this.Ka < SMALLVALUE ? "NA" : +this.Ka.toExponential(5))
		//Ks
		result.push(this.Ks < SMALLVALUE ? "NA" : +this.Ks.toExponential(5))
		//Ka/Ks
		result.push(this.Ks < SMALLVALUE || this.Ks == NA || this.Ka == NA ? "NA" : +(this.Ka / this.Ks).toExponential(5));
		//Fisher's test: p_value
		result.push(this.Sd < SMALLVALUE || this.Nd < SMALLVALUE || this.S < SMALLVALUE || this.N < SMALLVALUE ? "NA" : +this.fisher(this.Sd, this.Nd, this.S - this.Sd, this.N - this.Nd).toExponential(5));
		//Length of compared pairwise sequences
		result.push(length)
		//Synonymous(S) sites
		result.push(this.S < SMALLVALUE ? "NA" : +this.S.toExponential(5))
		//Nonsynonymous(N) sites
		result.push(this.N < SMALLVALUE ? "NA" : +this.N.toExponential(5))
		//L[0], L[2], L[4] only for Prof.Li's series(LWL85, LPB93...)
		result.push(this.L[0] < SMALLVALUE && this.L[2] < SMALLVALUE && this.L[4] < SMALLVALUE ? "NA" : [+this.L[0].toExponential(5), +this.L[2].toExponential(5), +this.L[4].toExponential(5)].join(":"))
		//Substitutions
		result.push(this.snp)
		//Sysnonymous Substitutions(Sd)
		result.push(this.Sd < SMALLVALUE ? "NA" : +this.Sd.toExponential(5))
		//Nonsysnonymous Substitutions(Nd)
		result.push(this.Nd < SMALLVALUE ? "NA" : +this.Nd.toExponential(5))

		//Si for Li's series' methods(LWL85, LPB93...)
		result.push(this.Si[0] < SMALLVALUE && this.Si[2] < SMALLVALUE && this.Si[4] < SMALLVALUE ? "NA" : [+this.Si[0].toExponential(5), +this.Si[2].toExponential(5), +this.Si[4].toExponential(5)].join(":"))

		//Vi for Li's series' methods(LWL85, LPB93...)
		result.push(this.Vi[0] < SMALLVALUE && this.Vi[2] < SMALLVALUE && this.Vi[4] < SMALLVALUE ? "NA" : [+this.Vi[0].toExponential(5), +this.Vi[2].toExponential(5), +this.Vi[4].toExponential(5)].join(":"))
		//Divergence time or distance t = (S*Ks+N*Ka)/(S+N)
		result.push(this.t < SMALLVALUE ? "NA" : +this.t.toExponential(5))
		//Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)
		result.push(this.KAPPA.map(x=>+x.toExponential(5)).join(":"))
		//GC Content
		result.push(+GC[0].toExponential(5) + "(" + GC.map(x=>+x.toExponential(5)).slice(1, 4).join(":") + ")")

		//Maximum Likelihood Value
		result.push(this.lnL == NA ? "NA" : +this.lnL.toExponential(5));
		//AICc
		result.push(this.AICc == NA ? "NA" : +this.AICc.toExponential(5));
		//Akaike weight in model selection
		result.push(this.AkaikeWeight == NA ? "NA" : this.AkaikeWeight);
		//Selected Model according to AICc
		result.push(this.model == "" ? "NA" : this.model);

		/*
		//Standard Errors
		result.push(this.SEKa == NA ? "NA" : this.SEKa);
		result.push(this.SEKs == NA ? "NA" : this.SEKs);
		*/
		return result.join("\t") + "\n";
	}

	/**************************************************
	* Function: fisher
	* Input Parameter: double, double, double, double
	* Output: Compute p-value by Fisher exact test
	* Return Value: double
	***************************************************/
	fisher(sd, nd, s, n) {
		let denominator, numerator, prob_total, prob_current, sum, fac_sum;
		let matrix = [0, 0, 0, 0], R = [0, 0], C = [0, 0];
		let i, j;

		denominator = numerator = prob_total = prob_current = sum = fac_sum = 0.0;

		matrix[0] = sd; matrix[1] = s; matrix[2] = nd; matrix[3] = n;
		//Row & Column
		R[0] = matrix[0] + matrix[2]; R[1] = matrix[1] + matrix[3];
		C[0] = matrix[0] + matrix[1]; C[1] = matrix[2] + matrix[3];
		sum = R[0] + R[1];

		//Calculate the numberator that is a constant
		numerator += this.factorial(R[0]);
		numerator += this.factorial(R[1]);
		numerator += this.factorial(C[0]);
		numerator += this.factorial(C[1]);

		//Log of Factorial of N
		fac_sum = this.factorial(sum);
		for (i = 0, denominator = fac_sum; i < 4; i++) {
			denominator += this.factorial(matrix[i]);
		}
		//Probability of current situtation
		prob_current = Math.exp(numerator - denominator);

		//Two-tail probabilities if less than prob_current
		for (i = 0; i < R[0] + SMALLVALUE; i++) {
			matrix[0] = i;
			matrix[1] = C[0] - i;
			matrix[2] = R[0] - i;
			matrix[3] = R[1] - C[0] + i;
			if (matrix[0] > SMALLVALUE && matrix[1] > SMALLVALUE && matrix[2] > SMALLVALUE && matrix[3] > SMALLVALUE) {
				for (j = 0, denominator = fac_sum; j < 4; j++) {
					denominator += this.factorial(matrix[j]);
				}
				let temp = numerator - denominator;
				temp = Math.exp(numerator - denominator);
				if (temp <= prob_current) {
					prob_total += temp;
				}
			}
		}

		return prob_total;
	}
	factorial(n) {
		let temp = 1.0;
		if (n > 0) {
			n = n + 1;
			let x = 0;
			x += 0.1659470187408462e-06 / (n + 7);
			x += 0.9934937113930748e-05 / (n + 6);
			x -= 0.1385710331296526 / (n + 5);
			x += 12.50734324009056 / (n + 4);
			x -= 176.6150291498386 / (n + 3);
			x += 771.3234287757674 / (n + 2);
			x -= 1259.139216722289 / (n + 1);
			x += 676.5203681218835 / (n);
			x += 0.9999999999995183;
			temp = Math.log(x) - 5.58106146679532777 - n + (n - 0.5) * Math.log(n + 6.5);
		}

		return (temp);
	}
}
