data {
  int<lower=0> nAllRows;				// all observations
	int<lower=0> nEvalRows;				// observations excluding first captures for CJS
	int<lower=0> nYears;				// total years in the data
	int<lower=0> nRivers;				// total observed rivers in the data
	int<lower=0> nFirstObsRows;			// number of first obs. (num fish)
	int<lower=0> nLastPossibleRows;		// = nFirstObsRows
	int<lower=0> kStates;				// nRivers + 1
	int<lower=0> longestObsChain;		// longest run of observations for any fish
	int<lower=0> evalRows[ nEvalRows ];	// indicators for CJS observations
	int<lower=0> season[ nAllRows ];	// season of observation
	int<lower=0> year[ nAllRows ];		// year of observation
	int<lower=0> firstObsRows[ nFirstObsRows ];			// indicators of first obs.
	int<lower=0,upper=1> availableDATA[ nAllRows ]; 	// (1-eta[i,t]) : 1 if fish hasn't emigrated
	int<lower=0,upper=1> encDATA[ nAllRows];			// observation of capture in a given sampling
	int<lower=0> riverDATA[ nAllRows ];					// river of fish, if captured
	int<lower=0> lastPossibleRows[ nLastPossibleRows ];	// last entry for each fish chain
}


transformed data {
	int<lower=0,upper=kStates> Y[nAllRows];	// transformed observation variable
	int<lower=0,upper=kStates> X[nAllRows]; // transformed state variable

	/*
		X[i] = k =>	z[i]=0 / zRiv[i]=1 (dead river), 	if k=1
					z[i]=1 / zRiv[i]=k, 				1 < k <= kStates

		Y[i] = k =>	o[i]=0 / riverDATA[i]=0 (NA), 		if k=1
					o[i]=1 / riverDATA[i]=(k-1), 		1 < k <= kStates
	 */
	
	for (i in 1:nAllRows) {
		if (encDATA[i]==0 && riverDATA[i]==0) Y[i] <- 1;
		for (j in 2:kStates) {
			if (encDATA[i]==1 && riverDATA[i]==j-1) {
				Y[i] <- j;
				X[i] <- j;	// note that X is only set where Y is not 1
			}
		}
	}
}

parameters {
	real   pBeta[ 4, nYears, (nRivers+1) ];		// recapture hyperparam
	real phiBeta[ 4, nYears, (nRivers+1) ];		// survival hyperparam
	real psiBeta[ 4, nRivers, nRivers ];		// movement hyperparam
}

model {

	/* transformed parameters (but not needed for reporting) */
	real p[ nAllRows, nRivers ];					// recapture parameters
	real phi[ nAllRows, nRivers ];					// survival parameters
	vector[nRivers] psi[ nAllRows, nRivers ];		// movement parameters

	/* simplified HMM parameters */
	/* Describe the following transformed model:
		X[i] ~ categorical( theta_trans[i, X[i-1]] )
		Y[i] ~ categorical( theta_emit[ i, X[i]  ] )
	 */
	vector[kStates] theta_trans[nAllRows, kStates];	// hidden state transition params
	vector[kStates] theta_emit[nAllRows, kStates];	// state-obs emission params

	real ePsi[ nAllRows, nRivers, nRivers ];		// numerator for psi
	real lpsi[ nAllRows, nRivers, nRivers ];
	real sumPsi [nAllRows, nRivers];				// denominator for psi

	int eval_idx;		// index for evalRows
	int gamma_idx;		// index for gamma
	real acc[kStates];	// accumulator for forward algorithm
	real gamma[longestObsChain+1, kStates];	 // gamma[i, k] = log Pr(Y[1:i], X_known[1:i], X[i]=k)

	/*** Model Priors ***/

	/* Recapture priors */   // TODO Vectorize this
  	for( s in 1:4 ){    
		for(y in 1:nYears){
			for( r in 1:(nRivers+1) ){
				pBeta[ s,y,r ] ~ normal( 0,1.225 );
			}
		}
	}

	/* survival priors */
	for( s in 1:4 ){    
		for(y in 1:nYears){  
			for( r in 1:(nRivers+1) ){  
				phiBeta[ s,y,r ] ~ normal( 0,1.225 );
			}
		}
	}

	/* Psi Priors */
	for( s in 1:4 ) {
		for(r in 1:(nRivers)){
			for(r2 in 1:(nRivers)){
				psiBeta[ s,r,r2 ]~ normal(0,2.25);
			}
		}
	}

	/*** Recapture (p) Model ***/
	for (i in 1:nEvalRows) {
		eval_idx <- evalRows[i]+1;
		for (j in 1:nRivers+1) {
			p[ eval_idx, j ] <- inv_logit(pBeta[ season[ eval_idx ], year[ eval_idx ], j]);
		}
	}

	/*** Survival (phi) Model ***/
	for (i in 1:nEvalRows) {
		eval_idx <- evalRows[i];
		for (j in 1:nRivers+1) {
			phi[ eval_idx, j ] <- inv_logit( phiBeta[ season[ eval_idx ], year[ eval_idx ], j]);
		}
	}

	/*** Movement (psi) Model ***/
	for (i in 1:nEvalRows) {
		eval_idx <- evalRows[i];
		for (r1 in 1: nRivers) { 
			for (r2 in 1:nRivers) {
				// normal priors on logit
				lpsi[eval_idx,r1,r2] <- psiBeta[season[eval_idx],r1,r2];
				ePsi[eval_idx,r1,r2] <- exp(lpsi[eval_idx,r1, r2]) *  (1-(r1==r2));
			}
			// constrain each set of psi's to sum to one. TODO: replace with cleaner simplex logic?
			sumPsi[eval_idx,r1] <- sum(ePsi[eval_idx,r1]);
			for (r2 in 1:nRivers) {
				psi[eval_idx,r1,r2] <- 
						( ePsi[eval_idx,r1, r2] / (1+sumPsi[eval_idx,r1]) ) * (1 - (r1==r2)) + 
                  		( 1 / (1+sumPsi[eval_idx,r1]) )  *  (r1==r2);

			}
		}
	}

	/* Simplified HMM Parameters */
	for (i in 1:nEvalRows) {
		eval_idx <- evalRows[i]+1;
		for (j in 1:kStates) {
			for (k in 1:kStates) {
				theta_emit[eval_idx,j,k] <- 0;
				theta_trans[eval_idx,j,k] <- 0;
			}
		}

		/* emission probability */
		theta_emit[eval_idx,1,1] <- 1;
		for (k in 2:5) theta_emit[eval_idx,k,1] <- 1 - (p[eval_idx,k-1] * availableDATA[eval_idx]);	
		for (k in 2:5) theta_emit[eval_idx,k,k] <- p[eval_idx,k-1] * availableDATA[eval_idx];

		/* transition probability */
		theta_trans[eval_idx,1,1] <- 1;
		for (k in 2:5) theta_trans[eval_idx,k,1] <- 1 - phi[eval_idx-1,k-1];
		for (j in 2:5) {
			for (k in 2:5) {
				theta_trans[eval_idx,j,k] <- phi[eval_idx-1,j-1] * psi[eval_idx-1,j-1,k-1];
			}
		}
	}

	/*** Likelihood ***/
	// Do a single pass over the whole data, accumulating sums as necessary

	// one run of the forward algo per chain
	for (chain in 1:nFirstObsRows) {
		// initialize gamma vals
		for (i in 1:longestObsChain+1) {
			for (j in 1:kStates) {
				gamma[i,j] <- 0;
			}
		}

		for (i in (firstObsRows[chain]+1):lastPossibleRows[chain]) {
			gamma_idx <- i - firstObsRows[chain] + 1;
			if (Y[i] != 1) {
				// current and last obs known
				if (Y[i-1] != 1) {
					gamma[ gamma_idx, X[i] ] <- gamma[gamma_idx-1, X[i-1]] + 
												log(theta_emit[i, X[i], Y[i]]) + 
												log(theta_trans[i, X[i-1], X[i]]);
				} 
				// current obs known, last obs. unknown
				else { 
					for (k in 1:kStates) {
						acc[k] <- gamma[gamma_idx-1, k] + 
									log(theta_emit[i, X[i], Y[i]]) + 
									log(theta_trans[i, k, X[i]]);
					}
					// marginalize over last X
					gamma[ gamma_idx, X[i] ] <- log_sum_exp(acc);
				}			
			} else {
				// current obs unkown, last obs known
				if (Y[i-1] != 1) {
					for (k in 1:kStates) {
						gamma[ gamma_idx, k ] <- gamma[gamma_idx-1, X[i-1]] +
												log(theta_emit[i, k, Y[i]]) +
												log(theta_trans[i, X[i-1], k]);
					}
				} 
				// current and last obs unknown
				else {
					for (k in 1:kStates) {
						for (j in 1:kStates) {
							acc[j] <- gamma[gamma_idx-1, j] +
										log(theta_emit[i, k, Y[i]]) +
										log(theta_trans[i, j, k]);
						}
						// marginalize over last X for each current val of X
						gamma[ gamma_idx, k ] <- log_sum_exp(acc);
					}
				}
			}
		}
		if (Y[lastPossibleRows[chain]] != 1) { // if very last obs of chain is known
			increment_log_prob(gamma[ gamma_idx, X[lastPossibleRows[chain]] ]);
		} else {	// else marginalize
			increment_log_prob(log_sum_exp(gamma[gamma_idx]));
		}
	}     
}
