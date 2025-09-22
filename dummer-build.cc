// Author: Martin C. Frith 2025
//         Patrick Styll   2025
// SPDX-License-Identifier: BSD-3-Clause

#include <fstream>
#include <iostream>

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

// Copied from HMMER 3.4:
#define OPT_symfrac 0.50
#define OPT_ere_aa  0.59
#define OPT_ere_nt  0.62
#define OPT_esigma  45.0

// Baum Welch parameters:
#define OPT_bw_maxiter 1000
#define OPT_bw_maxDiff 1e-6

// down-scale probabilities by this amount, to delay overflow:
const double scale1 = 1.0 / (1<<30) / (1<<30) / (1<<30) / (1<<30);
const double scale = scale1 * scale1 * scale1 * scale1;

// From "Dirichlet mixtures: a method for improved detection of weak
// but significant protein sequence homology"
// K Sjolander, K Karplus, M Brown, R Hughey, A Krogh, IS Mian, D Haussler

// There are 2 versions of this paper floating about the internet,
// with slightly different mixture coefficients!

// This is the version used by HMMER 3.4, and it's different from the
// version at the journal website.

const double blocks9[] = {
  0.178091,
  0.270671, 0.039848, 0.017576, 0.016415, 0.014268, 0.131916, 0.012391,
  0.022599, 0.020358, 0.030727, 0.015315, 0.048298, 0.053803, 0.020662,
  0.023612, 0.216147, 0.147226, 0.065438, 0.003758, 0.009621,

  0.056591,
  0.021465, 0.010300, 0.011741, 0.010883, 0.385651, 0.016416, 0.076196,
  0.035329, 0.013921, 0.093517, 0.022034, 0.028593, 0.013086, 0.023011,
  0.018866, 0.029156, 0.018153, 0.036100, 0.071770, 0.419641,

  0.0960191,
  0.561459, 0.045448, 0.438366, 0.764167, 0.087364, 0.259114, 0.214940,
  0.145928, 0.762204, 0.247320, 0.118662, 0.441564, 0.174822, 0.530840,
  0.465529, 0.583402, 0.445586, 0.227050, 0.029510, 0.121090,

  0.0781233,
  0.070143, 0.011140, 0.019479, 0.094657, 0.013162, 0.048038, 0.077000,
  0.032939, 0.576639, 0.072293, 0.028240, 0.080372, 0.037661, 0.185037,
  0.506783, 0.073732, 0.071587, 0.042532, 0.011254, 0.028723,

  0.0834977,
  0.041103, 0.014794, 0.005610, 0.010216, 0.153602, 0.007797, 0.007175,
  0.299635, 0.010849, 0.999446, 0.210189, 0.006127, 0.013021, 0.019798,
  0.014509, 0.012049, 0.035799, 0.180085, 0.012744, 0.026466,

  0.0904123,
  0.115607, 0.037381, 0.012414, 0.018179, 0.051778, 0.017255, 0.004911,
  0.796882, 0.017074, 0.285858, 0.075811, 0.014548, 0.015092, 0.011382,
  0.012696, 0.027535, 0.088333, 0.944340, 0.004373, 0.016741,

  0.114468,
  0.093461, 0.004737, 0.387252, 0.347841, 0.010822, 0.105877, 0.049776,
  0.014963, 0.094276, 0.027761, 0.010040, 0.187869, 0.050018, 0.110039,
  0.038668, 0.119471, 0.065802, 0.025430, 0.003215, 0.018742,

  0.0682132,
  0.452171, 0.114613, 0.062460, 0.115702, 0.284246, 0.140204, 0.100358,
  0.550230, 0.143995, 0.700649, 0.276580, 0.118569, 0.097470, 0.126673,
  0.143634, 0.278983, 0.358482, 0.661750, 0.061533, 0.199373,

  0.234585,
  0.005193, 0.004039, 0.006722, 0.006121, 0.003468, 0.016931, 0.003647,
  0.002184, 0.005019, 0.005990, 0.001473, 0.004158, 0.009055, 0.003630,
  0.006583, 0.003172, 0.003690, 0.002967, 0.002772, 0.002686,
};

// This Dirichlet mixture & description were copied from HMMER 3.4
// It doesn't satisfy Chargaff parity (A=T, C=G): undesirable for DNA?

// Match emission priors are trained on Rmark3 database
// Xref: ~wheelert/notebook/2011/0325_nhmmer_new_parameters

const double wheeler4[] = {
  0.24,
  0.16,  0.45,  0.12,   0.39,

  0.26,
  0.09,  0.03,  0.09,   0.04,

  0.08,
  1.29,  0.40,  6.58,   0.51,

  0.42,
  1.74,  1.49,  1.57,   1.95,
};

/* ----------from utils---------- */

bool isDash(const char *text) {
  return text[0] == '-' && text[1] == 0;
}

std::istream &fail(std::istream &s, const char *message) {
  std::cerr << message << "\n";
  s.setstate(std::ios::failbit);
  return s;
}

std::istream &openFile(std::ifstream &file, const char *name) {
  if (isDash(name)) return std::cin;
  file.open(name);
  if (!file) std::cerr << "can't open file: " << name << "\n";
  return file;
}

int badOpt() {
  std::cerr << "bad option value\n";
  return 1;
}

/* ------------------------------ */

/* ----------from priors--------- */

struct GapPriors {
  double match, insStart, delStart;
  double insEnd, insExtend;
  double delEnd, delExtend;
};

/* ------------------------------ */

struct DirichletMixture {
  const double *params;
  int componentCount;
};

struct MultipleAlignment {
  int sequenceCount;
  int alignmentLength;
  std::vector<unsigned char> alignment;
  std::string name;
};

void setCharToNumber(unsigned char *charToNumber, const char *alphabet) {
  for (int i = 0; alphabet[i]; ++i) {
    int c = alphabet[i];
    charToNumber[toupper(c)] = charToNumber[tolower(c)] = i;
  }
}

void normalize(double *x, int n) {
  double s = std::accumulate(x, x + n, 0.0);
  assert(s > 0);
  for (int i = 0; i < n; ++i) x[i] /= s;
}

double geometricMean(const double *values, int length, int step) {
  double s = 0;
  for (int i = 0; i < length; ++i) s += log(values[i * step]);
  return exp(s / length);
}

void setBackgroundProbs(double *bgProbs, const double *probs,
			int alphabetSize, int profileLength) {
  for (int i = 0; i < alphabetSize; ++i) {
    bgProbs[i] = geometricMean(probs + i, profileLength, 7 + alphabetSize);
  }
  normalize(bgProbs, alphabetSize);
}

// Probabilities are estimated from counts and a prior probability
// distribution (pseudocounts).

// If we want probabilities fitted to the input multiple alignment, we
// should use option --enone: that makes maxCountSum effectively
// infinite, so it has no effect.

// However, we often want probabilities fitted to more-distantly
// related sequences than those in the input.  To do that, we
// downscale the counts (relative to the prior) when the total count
// is high (> maxCountSum).  This is slightly different from
// (supposedly, better than) HMMER's "entropy weighting".

double getProb(double count1, double count2,
	       double pseudocount1, double pseudocount2, double maxCountSum) {
  double s = count1 + count2;
  double r = (s > maxCountSum) ? maxCountSum / s : 1.0;
  return (count1 * r + pseudocount1) / (s * r + (pseudocount1 + pseudocount2));
}

void applyDirichletMixture(DirichletMixture dmix,
			   int alphabetSize, double maxCountSum,
			   const double *counts, double *probEstimate) {
  int componentSize = 1 + alphabetSize;  // 1 mixture coefficient + alphas
  double countSum = std::accumulate(counts, counts + alphabetSize, 0.0);
  double rescale = (countSum > maxCountSum) ? maxCountSum / countSum : 1.0;
  std::vector<double> logs(dmix.componentCount);
  std::fill_n(probEstimate, alphabetSize, 0.0);

  for (int i = 0; i < dmix.componentCount; ++i) {
    const double *alphas = dmix.params + i * componentSize + 1;
    double alphaSum = std::accumulate(alphas, alphas + alphabetSize, 0.0);
    logs[i] = lgamma(alphaSum) - lgamma(countSum * rescale + alphaSum + 1);
    for (int j = 0; j < alphabetSize; ++j) {
      logs[i] += lgamma(counts[j] * rescale + alphas[j]) - lgamma(alphas[j]);
    }
  }

  double maxLog = *std::max_element(logs.begin(), logs.end());

  for (int i = 0; i < dmix.componentCount; ++i) {
    double componentProb = dmix.params[i * componentSize];
    const double *alphas = dmix.params + i * componentSize + 1;
    double mul = componentProb * exp(logs[i] - maxLog);
    for (int j = 0; j < alphabetSize; ++j) {
      probEstimate[j] += mul * (counts[j] * rescale + alphas[j]);
    }
  }

  normalize(probEstimate, alphabetSize);
}

void countsToParamCounts(const double* counts, double* paramCounts,
       int alphabetSize, int profileLength) {
  int countsPerPosition = 7 + alphabetSize;

  for (int i = 0; i <= profileLength; i++) {
    const double *c = counts + i * countsPerPosition;
    double *p = paramCounts + i * countsPerPosition;

    // copy letter counts
    std::copy(c + 7, c + 7 + alphabetSize, p + 7);

    // copy gap counts
    p[0] = c[1]; // gamma
    p[1] = c[2]; // delta
    p[2] = c[3]; // alpha
    p[3] = c[3]; // betap
    p[4] = c[4]; // beta
    p[5] = c[5]; // epsilonp
    p[6] = c[6]; // epsilon
  }
}

void gapCountsToProbs(const GapPriors &gp, double maxCountSum,
		      const double *counts, double *probs) {
  double alnBeg = counts[0]; // etap
  double match  = counts[1]; //gamma
  double delBeg = counts[2]; //delta
  double insBeg = counts[3]; //alpha
  double insEnd = insBeg;    //betap
  double insExt = counts[4]; //beta
  double delEnd = counts[5]; //epsilonp
  double delExt = counts[6]; //epsilon

  // The "- alnBeg" avoids the insertion start probability changing
  // if the aligned sequences are reversed:
  double notIns = match + delBeg - alnBeg;

  double gpNotIns = gp.match + gp.delStart;
  double a = getProb(insBeg, notIns, gp.insStart, gpNotIns, maxCountSum);
  probs[1] = a;  // insertion start probability

  double b = getProb(insExt, insEnd, gp.insExtend, gp.insEnd, maxCountSum);
  probs[3] = 1 - b;
  probs[4] = b;  // insertion extend probability

  double d = getProb(delBeg, match, gp.delStart, gp.match, maxCountSum);
  probs[0] = (1 - a) * (1 - d);
  probs[2] = (1 - a) * d;  // deletion start probability

  double e = getProb(delExt, delEnd, gp.delExtend, gp.delEnd, maxCountSum);
  probs[5] = 1 - e;
  probs[6] = e;  // deletion extend probability
}

void countsToProbs(DirichletMixture dmix, const GapPriors &gp,
		   int alphabetSize, double maxCountSum,
		   int profileLength, const double *counts, double *probs) {
  int countsPerPosition = 7 + alphabetSize;

  for (int i = 0; ; ++i) {
    const double *c = counts + i * countsPerPosition;
    /* */ double *p = probs  + i * countsPerPosition;

    gapCountsToProbs(gp, maxCountSum, c, p);

    if (i == profileLength) break;

    applyDirichletMixture(dmix, alphabetSize, maxCountSum, c + 7, p + 7);
  }

  // These gap probabilities are never used - set them similarly to HMMER:

  probs[5] = 1;  // delEnd
  probs[6] = 0;  // delExtend

  double *p = probs + profileLength * countsPerPosition;
  p[0] = 1 - p[1];  // match
  p[2] = 0;         // delStart

  setBackgroundProbs(p + 7, probs + 7, alphabetSize, profileLength);
}

std::istream &checkAlignmentBlock(std::istream &in, MultipleAlignment &ma,
				  int lineCount) {
  if (lineCount) {
    if (ma.sequenceCount && lineCount != ma.sequenceCount) {
      return fail(in, "alignment blocks should have same number of sequences");
    }
    ma.sequenceCount = lineCount;
  }
  return in;
}

std::istream &readMultipleAlignment(std::istream &in, MultipleAlignment &ma) {
  std::vector<std::string> lines;
  std::string line, word;
  int lineCount = 0;
  ma.sequenceCount = 0;
  ma.alignmentLength = 0;
  ma.name = "Anon";

  while (getline(in, line)) {
    std::istringstream iss(line);
    iss >> word;
    if (!iss) {
      if (!checkAlignmentBlock(in, ma, lineCount)) return in;
      lineCount = 0;
    } else if (word == "#=GF") {
      iss >> word;
      if (word == "ID") {
	iss >> ma.name;
      }
    } else if (word == "//") {
      break;
    } else if (word[0] != '#') {
      iss >> word;
      if (!iss) return fail(in, "bad alignment file");
      if (lineCount) {
	if (word.size() != lines.back().size()) {
	  return fail(in, "aligned sequences should have the same length");
	}
      } else {
	ma.alignmentLength += word.size();
      }
      lines.push_back(word);
      ++lineCount;
    }
  }

  if (!checkAlignmentBlock(in, ma, lineCount)) return in;

  ma.alignment.resize(ma.sequenceCount * ma.alignmentLength);
  unsigned char *a = ma.alignment.data();

  for (int i = 0; i < ma.sequenceCount; ++i) {
    for (size_t j = 0; j < lines.size(); j += ma.sequenceCount) {
      for (size_t k = 0; k < lines[i + j].size(); ++k) {
	unsigned char c = lines[i + j][k];
	*a++ = c;
      }
    }
  }

  return in;
}

bool isProteinAlignment(const MultipleAlignment &ma) {
  int counts[256] = {0};
  for (auto x : ma.alignment) ++counts[x];

  int tot = 0;
  for (int x = 'A'; x <= 'Z'; ++x) {
    int y = tolower(x);
    tot += counts[x] + counts[y];
  }

  int dna = 0;
  const char dnaChars[] = "ACGTNU";
  for (int i = 0; dnaChars[i]; ++i) {
    int x = dnaChars[i];
    int y = tolower(x);
    dna += counts[x] + counts[y];
  }

  return 1.0 * dna / tot < 0.9;  // xxx ???
}

void markEndGaps(MultipleAlignment &ma, int alphabetSize) {
  const int midGap = alphabetSize + 1;
  const int endGap = alphabetSize + 2;

  for (int i = 0; i < ma.sequenceCount; ++i) {
    unsigned char *seq = ma.alignment.data() + i * ma.alignmentLength;
    for (int j = 0; j < ma.alignmentLength && seq[j] == midGap; ++j) {
      seq[j] = endGap;
    }
    for (int j = ma.alignmentLength; j > 0 && seq[j-1] == midGap; --j) {
      seq[j-1] = endGap;
    }
  }
}

void makeSequenceWeights(const MultipleAlignment &ma, int alphabetSize,
			 double symfrac, double *weights) {
  const int midGap = alphabetSize + 1;
  const int endGap = alphabetSize + 2;
  std::vector<int> positionCounts(ma.sequenceCount);

  for (int i = 0; i < ma.alignmentLength; ++i) {
    int counts[32] = {0};
    const unsigned char *seq = ma.alignment.data() + i;
    for (int j = 0; j < ma.sequenceCount; ++j) {
      ++counts[seq[j * ma.alignmentLength]];
    }
    int totalCount = ma.sequenceCount - counts[endGap];
    int nonGapCount = totalCount - counts[midGap];

    if (nonGapCount > 0 && nonGapCount >= symfrac * totalCount) {
      int types = 0;
      for (int k = 0; k < alphabetSize; ++k) {
	types += (counts[k] > 0);
      }
      for (int j = 0; j < ma.sequenceCount; ++j) {
	int x = seq[j * ma.alignmentLength];
	if (x < alphabetSize) {
	  weights[j] += 1.0 / (types * counts[x]);
	  ++positionCounts[j];
	}
      }
    }
  }

  for (int j = 0; j < ma.sequenceCount; ++j) {
    if (positionCounts[j]) weights[j] /= positionCounts[j];
  }

  double maxWeight = *std::max_element(weights, weights + ma.sequenceCount);
  if (maxWeight) {
    for (int j = 0; j < ma.sequenceCount; ++j) weights[j] /= maxWeight;
  }
}

void countEvents(const MultipleAlignment &ma, int alphabetSize, double symfrac,
		 const double *weights, double weightSum,
		 std::vector<double> &allCounts, std::vector<int> &columns) {
  const int midGap = alphabetSize + 1;
  const int endGap = alphabetSize + 2;
  std::vector<char> states(ma.sequenceCount);
  double alnBeg = 0;
  double insBeg = 0;
  double insExt = 0;
  double delEnd = 0;

  for (int i = 0; i < ma.alignmentLength; ++i) {
    double counts[32] = {0};
    const unsigned char *seq = ma.alignment.data() + i;
    for (int j = 0; j < ma.sequenceCount; ++j) {
      counts[seq[j * ma.alignmentLength]] += weights[j];
    }
    double totalCount = weightSum - counts[endGap];
    double nonGapCount = totalCount - counts[midGap];

    if (nonGapCount > 0 && nonGapCount >= symfrac * totalCount) {
      columns.push_back(i);
      double delBeg = 0;
      double delExt = 0;
      for (int j = 0; j < ma.sequenceCount; ++j) {
	int x = seq[j * ma.alignmentLength];
	if (x <= alphabetSize) {  // an aligned letter in a "match" column
	  if (states[j] ==  0 )	alnBeg += weights[j];
	  if (states[j] == 'd') delEnd += weights[j];
	  states[j] = 'm';
	} else if (x == midGap) {  // a gap in a "match" column
	  if (states[j] != 'd') delBeg += weights[j];
	  // xxx assumes extension (not restart) of a deletion:
	  if (states[j] == 'd') delExt += weights[j];
	  states[j] = 'd';
	}
      }
      allCounts.push_back(alnBeg);
      allCounts.push_back(nonGapCount);
      allCounts.push_back(delBeg);
      allCounts.push_back(insBeg);  // insEnd = insBeg
      allCounts.push_back(insExt);
      allCounts.push_back(delEnd);
      allCounts.push_back(delExt);
      allCounts.insert(allCounts.end(), counts, counts + alphabetSize);
      alnBeg = insBeg = insExt = delEnd = 0;
    } else {  // this position is defined as "insertion"
      for (int j = 0; j < ma.sequenceCount; ++j) {
	int x = seq[j * ma.alignmentLength];
	if (x <= alphabetSize) {  // a non-gap symbol
	  if (states[j] ==  0 ) alnBeg += weights[j];
	  if (states[j] == 'd') delEnd += weights[j];
	  if (states[j] != 'i') insBeg += weights[j];
	  // xxx assumes extension (not restart) of an insertion:
	  if (states[j] == 'i') insExt += weights[j];
	  states[j] = 'i';
	}
      }
    }
  }

  allCounts.insert(allCounts.end(), 7, 0.0);  // final gap counts must all be 0

  // We'll estimate probabilities from these counts, but there are 2 problems:
  // 1. Adjacent inserts (or deletes) are assumed to be extension, not restart.
  // 2. The counts assume the alignment is exactly correct and certain.
  // We will mitigate both problems by using a Baum-Welch algorithm.
}

void getSequenceWithoutGaps(
     const MultipleAlignment &ma, int alphabetSize,
     std::vector<unsigned char> &seqNoGap, std::vector<int> &seqLengths) {
  for (int i = 0; i < ma.sequenceCount; i++) {
    const unsigned char *seq = ma.alignment.data() + i * ma.alignmentLength;
    for (int j = 0; j < ma.alignmentLength; j++) {
      if (seq[j] <= alphabetSize) {
        seqNoGap.push_back(seq[j]);
        seqLengths[i]++;
      }
    }
  }
  seqNoGap.push_back(0);  // so it's safe to read one letter past the end
}

void forward(
     unsigned char *seq, int seqLength, std::vector<double> &probs,
     int profileLength, int width, double *wSum,
     std::vector<double> &X, std::vector<double> &Y, std::vector<double> &Z) {

  int cols = seqLength + 2;
  double aPrime, bPrime, dPrime, ePrime;  // parameters for BW
  double alphaProb, betaProb, deltaProb, epsilonProb, epsilonProb1;

  std::fill_n(&X[0], seqLength+1, 0.0);
  std::fill_n(&Y[1], seqLength+1, 0.0);

  for (int i = 1; i <= profileLength+1; i++) {

    alphaProb    = probs[(i-1) * width + 1];  // insStart
    betaProb     = probs[(i-1) * width + 4];  // insExt
    deltaProb    = probs[(i-1) * width + 2];  // delStart
    epsilonProb  = probs[(i-1) * width + 6];  // delExt

    epsilonProb1 = (i == profileLength+1) ? 0.0 : probs[i * width + 6];

    aPrime = alphaProb * (1.0 - betaProb);
    bPrime = betaProb;

    // can use arbitrary values for S_m, d_m, e_m
    dPrime = deltaProb   * (1 - epsilonProb1);
    ePrime = epsilonProb * (1 - epsilonProb1) / (1 - epsilonProb);
    if (!std::isfinite(ePrime)) ePrime = 0.0; // avoid NaN issues

    X[(i-1) * cols] = 0.0;
    double z = Z[i * cols] = 0.0;

    for (int j = 1; j <= seqLength+1; j++) {

      int letter = seq[j-1];

      // letter probs / background probs
      double S = (1 - alphaProb - deltaProb)
               * probs[(i-1) * width + (7 + letter)]
               / probs[profileLength * width + (7 + letter)];
      if (!std::isfinite(S)) S = 0.0; // if letter has 0 background probability

      double y = Y[(i-1) * cols + j];
      double w = X[(i-1) * cols + (j-1)] + y + z + scale;
      z = aPrime * w + bPrime * z;
      *wSum += w;

      X[i * cols + j] = S * w;
      Y[i * cols + j] = dPrime * w + ePrime * y;
      Z[i * cols + j] = z;
    }
  }
}

void backward(unsigned char *seq, int seqLength,
     std::vector<double> &probs, int profileLength, int width,
     std::vector<double> &Wbar, std::vector<double> &Ybar, std::vector<double> &Zbar) {

  int cols = seqLength + 2;
  double aPrime, bPrime, dPrime, ePrime;  // parameters for BW
  double alphaProb, betaProb, deltaProb, epsilonProb, epsilonProb1;

  std::fill_n(&Wbar[(profileLength+1) * cols + 1], seqLength+1, 0.0);
  std::fill_n(&Ybar[(profileLength+1) * cols + 0], seqLength+1, 0.0);

  for (int i = profileLength; i >= 0; i--) {

    alphaProb    = probs[i * width + 1];  // insStart
    betaProb     = probs[i * width + 4];  // insExt
    deltaProb    = probs[i * width + 2];  // delStart
    epsilonProb  = probs[i * width + 6];  // delExt

    epsilonProb1 = (i == profileLength) ? 0.0 : probs[(i+1) * width + 6];

    aPrime = alphaProb * (1.0 - betaProb);
    bPrime = betaProb;

    // can use arbitrary values for S_m, d_m, e_m
    dPrime = deltaProb   * (1 - epsilonProb1);
    ePrime = epsilonProb * (1 - epsilonProb1) / (1 - epsilonProb);
    if (!std::isfinite(ePrime)) ePrime = 0.0; // avoid NaN issues

    Wbar[(i+2) * cols - 1] = 0.0;
    double z = Zbar[(i+1) * cols - 1] = 0.0;

    for (int j = seqLength; j >= 0; j--) {

      int letter = seq[j];

      // letter probs / background probs
      double S = (1 - alphaProb - deltaProb)
               * probs[i * width + (7 + letter)]
               / probs[profileLength * width + (7 + letter)];
      if (!std::isfinite(S)) S = 0.0; // if letter has 0 background probability

      double x = S * Wbar[(i+1) * cols + (j+1)];
      double y = Ybar[(i+1) * cols + j];
      double w = x + dPrime * y + aPrime * z + scale;
      z = w + bPrime * z;

      Wbar[i * cols + j] = w;
      Ybar[i * cols + j] = w + ePrime * y;
      Zbar[i * cols + j] = z;
    }
  }
}

void calculateTransitionCounts(
     std::vector<double> &counts, int profileLength, int seqLength,
     std::vector<double> &Wbar, std::vector<double> &Ybar,
     std::vector<double> &Zbar, std::vector<double> &X,
     std::vector<double> &Y,    std::vector<double> &Z,
     double wt, int width, int alphabetSize, unsigned char *seq) {

  int cols = seqLength + 2;
  double alpha, beta, delta, epsilon;
  double gamma, betap, epsilonp, etap;

  // estimate counts for every position
  for (int i = 1; i <= profileLength+1; i++) {

    double emis[32] = {0};
    betap    = 0.0;  // expected count of (1 - beta)
    beta     = 0.0;
    epsilonp = 0.0;  // expected count of (1 - epsilon)
    epsilon  = 0.0;
    etap     = 0.0;  // expected count of (1 - eta)

    for (int j = 1; j <= seqLength+1; j++) {
      betap    += Z[i * cols + (j-1)] * Wbar[(i-1) * cols + (j-1)];
      beta     += Z[i * cols + (j-1)] * Zbar[(i-1) * cols + (j-1)];
      epsilonp += Y[(i-1) * cols + j] * Wbar[(i-1) * cols + (j-1)];
      epsilon  += Y[(i-1) * cols + j] * Ybar[(i-1) * cols + (j-1)];
      etap     += Wbar[(i-1) * cols + (j-1)];
      int letter = seq[j-1];
      emis[letter] += X[i * cols + j] * Wbar[i * cols + j];
    }

    beta  -= betap;
    if (beta < -__DBL_EPSILON__)
      std::cerr << "Warning: Negative beta count at position " << i << "\n";

    epsilon  -= epsilonp;
    if (epsilon < -__DBL_EPSILON__)
      std::cerr << "Warning: Negative epsilon count at position " << i << "\n";

    // expected count of alpha
    alpha = betap;

    gamma = std::accumulate(emis, emis + alphabetSize + 1, 0.0);

    // update the HMM parameters
    counts[(i-1) * width + 0] += etap * scale * wt;
    counts[(i-1) * width + 1] += gamma        * wt;

    counts[(i-1) * width + 3] += alpha        * wt;
    counts[(i-1) * width + 4] += beta         * wt;
    counts[(i-1) * width + 5] += epsilonp     * wt;
    counts[(i-1) * width + 6] += epsilon      * wt;

    // emission probabilities
    if (i == profileLength+1) continue; // no emissions from the end state
    for (int letter = 0; letter < alphabetSize; letter++) {
      counts[(i-1) * width + (7 + letter)] += emis[letter] * wt;
    }
  }

  for (int i = 1; i <= profileLength; i++) {
    delta = counts[i * width + 5] + counts[i * width + 6]
      - counts[(i-1) * width + 6];
    if (delta < -__DBL_EPSILON__)
      std::cerr << "Warning: Negative delta count at position " << i << "\n";
    counts[(i-1) * width + 2] = delta;
  }
}

int determineTerminationCondition(double bwMaxDiff,
       std::vector<double> &probsOld, std::vector<double> &probsNew) {
  double maxDiff = 0.0;
  for (size_t i = 0; i < probsOld.size(); i++) {
    double diff = fabs(probsOld[i] - probsNew[i]);
    if (diff > maxDiff) maxDiff = diff;
    if (std::isnan(probsNew[i])) return -1; // overflow
  }
  return (maxDiff <= bwMaxDiff) ? 1 : 0; // converged or not
}

void baumWelch(std::vector<double> &counts, const MultipleAlignment &ma,
     int alphabetSize, DirichletMixture dmix, const GapPriors &gp,
     int profileLength, const double *weights,
     double maxIter, double bwMaxDiff) {

  // transitions + alphabet
  int width = 7 + alphabetSize;
  std::vector<double> probsOld(counts.size() + alphabetSize);
  std::vector<double> probsNew(counts.size() + alphabetSize);

  int rows = profileLength      + 2;
  int cols = ma.alignmentLength + 2;

  std::vector<double> Wbar(rows * cols, 0.0);
  std::vector<double> Ybar(rows * cols, 0.0);
  std::vector<double> Zbar(rows * cols, 0.0);
  std::vector<double> X(rows * cols, 0.0);
  std::vector<double> Y(rows * cols, 0.0);
  std::vector<double> Z(rows * cols, 0.0);

  /* Get gapless sequences and their lengths. */
  std::vector<unsigned char> seqsNoGap;
  std::vector<int>  seqsLengths(ma.sequenceCount, 0);
  getSequenceWithoutGaps(ma, alphabetSize, seqsNoGap, seqsLengths);

  countsToProbs(dmix, gp, alphabetSize, 1e37, profileLength,
      counts.data(), probsOld.data());

  int termCond;
  int num_iterations = 0;

  do {

    /* We aggregate the expected counts over all sequences. */
    for (size_t i = 0; i < counts.size(); i++) counts[i] = 0.0;
    std::fill(probsNew.begin(), probsNew.end(), 0.0);

    unsigned char* seqNoGap = seqsNoGap.data();

    for (int idx = 0; idx < ma.sequenceCount; idx++) {
      int seqLength = seqsLengths[idx];

      /* Forward pass, calculate X, Y, Z and
         the aggregated v (i.e. sum of w-values). */
      double v = 0.0; // accumulate the sum of weights
      forward(seqNoGap, seqLength, probsOld,
        profileLength, width, &v, X, Y, Z);

      // numbers overflowed to infinity
      if (v > DBL_MAX || std::isnan(v)) break;

      /* Backward pass, calculate Wbar, Ybar, Zbar. */
      backward(seqNoGap, seqLength, probsOld,
        profileLength, width, Wbar, Ybar, Zbar);

      double wt = weights[idx] / (v * scale);

      /* Calculate and update parameters in new parameter counts. */
      calculateTransitionCounts(counts, profileLength, seqLength,
          Wbar, Ybar, Zbar, X, Y, Z, wt, width, alphabetSize, seqNoGap);

      seqNoGap += seqLength;
    }

    /* Turn the parameter counts into probabilities via priors. */
    countsToProbs(dmix, gp, alphabetSize, 1e37, profileLength,
        counts.data(), probsNew.data());

    /* Determine termination condition.
       Then overwrite the old probabilities with the new ones. */
    termCond = determineTerminationCondition(bwMaxDiff, probsOld, probsNew);

    if (termCond == -1) break; // overflow, keep old probs

    std::copy(probsNew.begin(), probsNew.end(), probsOld.begin());

    num_iterations++;

  } while (!termCond && num_iterations < maxIter);

}

void printProb(double prob) {
  if (prob > 0) {
    std::cout << "  " << abs(log(prob));
  } else {
    std::cout << "        *";
  }
}

std::vector<std::vector<double>> build_hmm(const char* filename, double symfrac,
          double ere, double esigma,
          double bwMaxiter, double bwMaxDiff,
          bool countOnly, GapPriors gp) {

  std::ifstream file;
  std::istream &in = openFile(file, filename);
  if (!file) exit(1);

  unsigned char charToNumber[256];
  MultipleAlignment ma;

  std::cout << std::fixed;

  std::vector<std::vector<double>> all_counts;

  while (readMultipleAlignment(in, ma)) {
    bool isProtein = isProteinAlignment(ma);
    const char *alphabet = isProtein ? "ACDEFGHIKLMNPQRSTVWY" : "ACGT";
    int alphabetSize = strlen(alphabet);
    memset(charToNumber, alphabetSize, 256);
    setCharToNumber(charToNumber, alphabet);
    if (!isProtein) setCharToNumber(charToNumber, "ACGU");
    charToNumber['-'] = charToNumber['.']
      = charToNumber['_'] = charToNumber['~'] = alphabetSize + 1;

    int width = 7 + alphabetSize;

    DirichletMixture dmix;
    dmix.params = isProtein ? blocks9 : wheeler4;
    dmix.componentCount = (isProtein ? sizeof blocks9 : sizeof wheeler4)
      / sizeof(double) / (1 + alphabetSize);

    for (auto &i : ma.alignment) i = charToNumber[i];

    markEndGaps(ma, alphabetSize);

    std::vector<double> weights(ma.sequenceCount);
    makeSequenceWeights(ma, alphabetSize, symfrac, weights.data());
    double weightSum = std::accumulate(weights.begin(), weights.end(), 0.0);

    std::vector<double> counts;
    std::vector<int> columns;
    countEvents(ma, alphabetSize, symfrac, weights.data(), weightSum,
        counts, columns);
    int profileLength = columns.size();

    if (!countOnly) {
      baumWelch(counts, ma, alphabetSize, dmix, gp,
        profileLength, weights.data(), bwMaxiter, bwMaxDiff);
    }

    /* As output, we will have the transition counts only. */
    std::vector<double> paramCounts(counts.size() + alphabetSize);
    countsToParamCounts(counts.data(), paramCounts.data(),
      alphabetSize, profileLength);

    for (int i = 0; i <= profileLength; i++) {
      std::vector<double> transitions_i(7, 0.0);
      for (int t = 0; t < 7; t++)
        transitions_i[t] = paramCounts[i * width + t];
      all_counts.push_back(std::move(transitions_i));
    }

  }

  return all_counts;
}

double relativeEntropy(const double *probs,
		       int alphabetSize, int profileLength) {
  int probsPerPosition = 7 + alphabetSize;
  const double *bgProbs = probs + profileLength * probsPerPosition;

  double r = 0;
  for (int i = 0; i < profileLength; ++i) {
    for (int j = 0; j < alphabetSize; ++j) {
      r += probs[j] * log2(probs[j] / bgProbs[j]);
    }
    probs += probsPerPosition;
  }
  return r;
}

std::vector<double> get_relative_entropies(const char* filename, double symfrac,
          double ere,
          double bwMaxiter, double bwMaxDiff,
          bool countOnly, GapPriors gp) {

  std::ifstream file;
  std::istream &in = openFile(file, filename);
  if (!file) exit(1);

  unsigned char charToNumber[256];
  MultipleAlignment ma;

  std::cout << std::fixed;

  std::vector<double> all_entropies;

  while (readMultipleAlignment(in, ma)) {
    bool isProtein = isProteinAlignment(ma);
    const char *alphabet = isProtein ? "ACDEFGHIKLMNPQRSTVWY" : "ACGT";
    int alphabetSize = strlen(alphabet);
    memset(charToNumber, alphabetSize, 256);
    setCharToNumber(charToNumber, alphabet);
    if (!isProtein) setCharToNumber(charToNumber, "ACGU");
    charToNumber['-'] = charToNumber['.']
      = charToNumber['_'] = charToNumber['~'] = alphabetSize + 1;

    DirichletMixture dmix;
    dmix.params = isProtein ? blocks9 : wheeler4;
    dmix.componentCount = (isProtein ? sizeof blocks9 : sizeof wheeler4)
      / sizeof(double) / (1 + alphabetSize);

    for (auto &i : ma.alignment) i = charToNumber[i];

    markEndGaps(ma, alphabetSize);

    std::vector<double> weights(ma.sequenceCount);
    makeSequenceWeights(ma, alphabetSize, symfrac, weights.data());
    double weightSum = std::accumulate(weights.begin(), weights.end(), 0.0);

    std::vector<double> counts;
    std::vector<int> columns;
    countEvents(ma, alphabetSize, symfrac, weights.data(), weightSum,
        counts, columns);
    int profileLength = columns.size();

    if (!countOnly) {
      baumWelch(counts, ma, alphabetSize, dmix, gp,
        profileLength, weights.data(), bwMaxiter, bwMaxDiff);
    }

    std::vector<double> probs(counts.size() + alphabetSize);
    countsToProbs(dmix, gp, alphabetSize, 1e37, profileLength,
        counts.data(), probs.data());

    /* As output, we will now have the relative entropies only (one per MSA). */
    double re = relativeEntropy(probs.data() + 7, alphabetSize, profileLength);

    all_entropies.push_back(re);
  }

  return all_entropies;
}

/* C-Types Wrapper */

extern "C" {

/* Need a struct to mimic Python lists */
struct ArrayDouble {
    double* data;
    size_t size;
};

ArrayDouble build_hmm_c(const char* filename,
                      double symfrac,
                      double ere, double esigma,
                      double bwMaxiter, double bwMaxDiff,
                      bool countOnly,
                      GapPriors gp)
{
    std::vector<std::vector<double>> probs =
        build_hmm(filename, symfrac, ere, esigma,
                  bwMaxiter, bwMaxDiff, countOnly, gp);
    
    /* Flatten probs into 1d array for Python... */

    size_t total_size = probs.size() * 7;
    double* res = static_cast<double*>(std::malloc(total_size * sizeof(double)));

    size_t idx = 0;
    for (const auto& v : probs) {
      std::copy(v.begin(), v.end(), res + idx);
      idx += v.size();
    }

    ArrayDouble arr;
    arr.data = res;
    arr.size = total_size;
    return arr;
}

ArrayDouble get_relative_entropies_c(const char* filename,
                      double symfrac,
                      double ere,
                      double bwMaxiter, double bwMaxDiff,
                      bool countOnly,
                      GapPriors gp)
{
    std::vector<double> entropies =
        get_relative_entropies(filename, symfrac, ere,
                  bwMaxiter, bwMaxDiff, countOnly, gp);
    
    size_t total_size = entropies.size();
    double* res = static_cast<double*>(std::malloc(total_size * sizeof(double)));
    std::copy(entropies.begin(), entropies.end(), res);

    ArrayDouble arr;
    arr.data = res;
    arr.size = total_size;
    return arr;

}

} // extern "C"

// This really only exists for testing purposes...
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <alignment_file>\n";
        return 1;
    }

    // Use default parameters for testing
    double symfrac   = OPT_symfrac;
    double ere       = OPT_ere_aa;
    double esigma    = OPT_esigma;
    double bwMaxiter = OPT_bw_maxiter;
    double bwMaxDiff = OPT_bw_maxDiff;
    bool countOnly = false;
    GapPriors gp = {0.9, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

    auto all_probs = build_hmm(argv[1], symfrac, ere, esigma, bwMaxiter, bwMaxDiff, countOnly, gp);

    std::cout << "Number of alignments processed: " << all_probs.size() << "\n";
    for (size_t i = 0; i < all_probs.size(); ++i) {
      std::cout << "Alignment " << i+1 << " mean transitions: ";
      for (size_t j = 0; j < all_probs[i].size(); ++j) {
        std::cout << all_probs[i][j] << " ";
      }
      std::cout << "\n";
    }
    return 0;
}