
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "LandauVishkin.h"

int InitiateLVCompute(short L[32][MAX_K+1][2*MAX_K+1], short thread_n)
{   
    int i,j, p_i;
	for(p_i = 0; p_i < thread_n; p_i++)
		for ( i = 0; i < MAX_K+1; i++) {
			for ( j = 0; j < 2*MAX_K+1; j++) {
				L[p_i][i][j] = -2;
			}
		}
    return 0;
}

int computeEditDistance(
        const char* text, int textLen, const char* pattern, int patternLen, int k, short L[MAX_K+1][2*MAX_K+1])
{
    
    _ASSERT(k < MAX_K);
    k = __min(MAX_K - 1, k); // enforce limit even in non-debug builds
    if (NULL == text) {
        // This happens when we're trying to read past the end of the genome.
        return -1;
    }

    const char* p = pattern;
    const char* t = text;
    int endl = __min(patternLen, textLen);
    const char* pend = pattern + endl; 
    while (p < pend) {
        
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, endl);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = endl;
done1:
    if (L[0][MAX_K] == endl) {
        int result = (patternLen > endl ? patternLen - endl : 0); // Could need some deletions at the end
        return result;
    }
    int e, d;
    for ( e = 1; e <= k; e++) {
        // Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
        for ( d = 0; d != e+1; d = (d > 0 ? -d : -d+1)) {
            int best = L[e-1][MAX_K+d] + 1; // up
            int left = L[e-1][MAX_K+d-1];
            if (left > best)
                best = left;
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best)
                best = right;

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t) { 
                int endl = __min(patternLen, textLen - d);
                const char* pend = pattern + endl;

                while (1) {
    
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, endl);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = endl;
                        break;
                    }
                    t += 8;
                }
            }

            if (best == patternLen) {
                return e;
            }
            L[e][MAX_K+d] = best;
        }
    }
    return -1;
}
bool writeCigar(char** o_buf, int* o_buflen, int count, char code)
{
    if (count <= 0) {
        return true;
    }
	
	/*
    if (format == EXPANDED_CIGAR_STRING) {
        int n = __min(*o_buflen, count);
        int i;
        for (i = 0; i < n; i++) {
            *(*o_buf)++ = code;
        }
        *o_buflen -= n;
        if (*o_buflen == 0) {
            *(*o_buf - 1) = '\0'; 
        }
        return *o_buflen > 0;
    } else if (format == COMPACT_CIGAR_STRING) {
        if (*o_buflen == 0) {
            *(*o_buf - 1) = '\0';
            return false;
        }
        int written = snprintf(*o_buf, *o_buflen, "%d%c", count, code);
        if (written > *o_buflen - 1) {
            *o_buf = '\0';
            return false;
        } else {
            *o_buf += written;
            *o_buflen -= written;
            return true;
        }
    } else if (format == COMPACT_CIGAR_BINARY) {
        // binary format with non-zero count byte followed by char (easier to examine programmatically)
        while (true) {
            if (*o_buflen < 3) {
                *(*o_buf) = '\0';
                return false;
            }
            *(*o_buf)++ = __min(count, 255);
            *(*o_buf)++ = code;
            *o_buflen -= 2;
            if (count <= 255) {
                return true;
            }
            count -= 255;
        }
    } else {
        printf("invalid cigar format %d\n", format);
        exit(1);
    }
	*/
	
	if (*o_buflen == 0) {
        *(*o_buf - 1) = '\0';
        return false;
    }
    int written = snprintf(*o_buf, *o_buflen, "%d%c", count, code);

	if (written > *o_buflen - 1) {
        *o_buf = '\0';
        return false;
    } else {
        *o_buf += written;
        *o_buflen -= written;
        return true;
    }		
}


int computeEditDistanceWithCigar(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, bool useM, 
    short L[MAX_K+1][2*MAX_K+1])//CigarFormat format, 
{
    _ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.

#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
	
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
        return 0;
    }
    int e;
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            int best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                int curD = d;
                int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
                    }
                    curE++;
                }
				if (useM && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
    return -1;
}

int computeEditDistance_s(
        const char* text, int textLen, const char* pattern, int patternLen, int k, short L[MAX_K+1][2*MAX_K+1], int* s_position)
{
    
    _ASSERT(k < MAX_K);
    k = __min(MAX_K - 1, k); // enforce limit even in non-debug builds
    if (NULL == text) {
        // This happens when we're trying to read past the end of the genome.
        return -1;
    }

    const char* p = pattern;
    const char* t = text;
    int endl = __min(patternLen, textLen);
    const char* pend = pattern + endl; 
    while (p < pend) {
        
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, endl);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = endl;
done1:
    if (L[0][MAX_K] == endl) {
        int result = (patternLen > endl ? patternLen - endl : 0); // Could need some deletions at the end
        return result;
    }

    int e, d;
	int best;
	int best_i;
    for ( e = 1; e <= k; e++) {
        // Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
        for ( d = 0; d != e+1; d = (d > 0 ? -d : -d+1)) {
            best = L[e-1][MAX_K+d] + 1; // up
            int left = L[e-1][MAX_K+d-1];
            if (left > best)
				best = left;
                
            int right = L[e-1][MAX_K+d+1] + 1;
			
            if (right > best)
				best = right;

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t) {
                int endl = __min(patternLen, textLen - d);
                const char* pend = pattern + endl;

                while (1) {
    
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, endl);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = endl;
                        break;
                    }
                    t += 8;
                }
            }

            if (best == patternLen) {
                return e;
            }
            L[e][MAX_K+d] = best;
        }
    }
	
	int max_L = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
			max_L = L[e - 1][MAX_K + best_i];
			
	(*s_position) = max_L;
	
    return -1;
}

int computeEditDistanceWithCigar_s(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1])//bool useM, CigarFormat format, 
{
    _ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
	
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					//printf("%s\n", cigarBuf);
					
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							
							//printf("%s\n", cigarBuf);
							
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}
							accumulatedMs = 0;
						}
						
						//printf("%s\n", cigarBuf);
						
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					
					//printf("%s\n", cigarBuf);
					
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	max_L = curD;
	
	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];
        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	accumulatedMs = L[0][MAX_K+0];
	
    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				
				//printf("%s\n", cigarBuf);
				
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}
				accumulatedMs = 0;
			}
			
			//printf("%s\n", cigarBuf);
			
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		
		//printf("%s\n", cigarBuf);
		
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - L[e - 1][MAX_K + max_L], 'S')) {//, format
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string

	
    //return -1;
	if(max_L < 0)	max_L = -max_L;
	
	return	max_L;
}

int computeEditDistanceWithCigar_s_nm(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], int* nm_score)//bool useM, CigarFormat format, 
{
    _ASSERT(k < MAX_K);
	
	int m_num = 0;
	
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
    char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}else{
								m_num += accumulatedMs;
							}
							accumulatedMs = 0;
						}
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}else{
						m_num += accumulatedMs;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
				
				(*nm_score) = patternLen - m_num;
				
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	int curD_re = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	curD_re = curD;
	//max_L = curD;
	
	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];
        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	accumulatedMs = L[0][MAX_K+0];
	
    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}else{
					m_num += accumulatedMs;
				}
				accumulatedMs = 0;
			}
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}else{
			m_num += accumulatedMs;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - max_L, 'S')) {//, format L[e - 1][MAX_K + max_L]
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string

	(*nm_score) = max_L - m_num;
	
    //return -1;
	if(curD_re < 0)	curD_re = -curD_re;
	
	return curD_re;
}

//for quality
int computeEditDistance_mis(
        const char* text, int textLen, const char* pattern, int patternLen, int k, short L[MAX_K+1][2*MAX_K+1], uint8_t* quality)
{
    _ASSERT(k < MAX_K);
    k = __min(MAX_K - 1, k); // enforce limit even in non-debug builds
    if (NULL == text) {
        // This happens when we're trying to read past the end of the genome.
        return -1;
    }
	
	uint8_t* quality_p = quality;
	
    const char* p = pattern;
    const char* t = text;
    int endl = __min(patternLen, textLen);
    const char* pend = pattern + endl; 
    while (p < pend) {
        
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
		x &= (*((_uint64*)quality_p));

        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, endl);
            goto done1;
        }
        p += 8;
        t += 8;
		quality_p += 8;
    }
    L[0][MAX_K] = endl;
done1:
    if (L[0][MAX_K] == endl) {
        int result = (patternLen > endl ? patternLen - endl : 0); // Could need some deletions at the end
        return result;
    }
    int e, d;
    for ( e = 1; e <= k; e++) {
        // Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
        for ( d = 0; d != e+1; d = (d > 0 ? -d : -d+1)) {
            int best = L[e-1][MAX_K+d] + 1; // up
            int left = L[e-1][MAX_K+d-1];
            if (left > best)
                best = left;
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best)
                best = right;

            const char* p = pattern + best;
            const char* t = (text + d) + best;
			quality_p = quality + best;
			
            if (*p == *t) { 
                int endl = __min(patternLen, textLen - d);
                const char* pend = pattern + endl;

                while (1) {
    
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
					x &= (*((_uint64*)quality_p));
					
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, endl);
                        break;
                    }
                    p += 8;
					quality_p += 8;
					
                    if (p >= pend) {
                        best = endl;
                        break;
                    }
                    t += 8;
                }
            }

            if (best == patternLen) {
                return e;
            }
            L[e][MAX_K+d] = best;
        }
    }
    return -1;
}




int computeEditDistance_mis_s(
        const char* text, int textLen, const char* pattern, int patternLen, int k, short L[MAX_K+1][2*MAX_K+1], uint8_t* quality, int16_t* s_position)
{
    _ASSERT(k < MAX_K);
    k = __min(MAX_K - 1, k); // enforce limit even in non-debug builds
    if (NULL == text) {
        // This happens when we're trying to read past the end of the genome.
        return -1;
    }
	
	uint8_t* quality_p = quality;
	
    const char* p = pattern;
    const char* t = text;
    int endl = __min(patternLen, textLen);
    const char* pend = pattern + endl; 
    while (p < pend) {
        
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
		x &= (*((_uint64*)quality_p));

        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, endl);
            goto done1;
        }
        p += 8;
        t += 8;
		quality_p += 8;
    }
    L[0][MAX_K] = endl;
done1:
    if (L[0][MAX_K] == endl) {
        int result = (patternLen > endl ? patternLen - endl : 0); // Could need some deletions at the end
        return result;
    }
    int e, d;
    for ( e = 1; e <= k; e++) {
        // Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
        for ( d = 0; d != e+1; d = (d > 0 ? -d : -d+1)) {
            int best = L[e-1][MAX_K+d] + 1; // up
            int left = L[e-1][MAX_K+d-1];
            if (left > best)
                best = left;
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best)
                best = right;

            const char* p = pattern + best;
            const char* t = (text + d) + best;
			quality_p = quality + best;
			
            if (*p == *t) { 
                int endl = __min(patternLen, textLen - d);
                const char* pend = pattern + endl;

                while (1) {
    
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
					x &= (*((_uint64*)quality_p));
					
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, endl);
                        break;
                    }
                    p += 8;
					quality_p += 8;
					
                    if (p >= pend) {
                        best = endl;
                        break;
                    }
                    t += 8;
                }
            }

            if (best == patternLen) {
                return e;
            }
            L[e][MAX_K+d] = best;
        }
    }
	
	int max_L = 0;
	int best_i = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
			max_L = L[e - 1][MAX_K + best_i];

	
	(*s_position) = patternLen - max_L;
	
	
    return -1;
}



int computeEditDistance_misboth(
        const char* text, int textLen, const char* pattern, int patternLen, int k, short L[MAX_K+1][2*MAX_K+1], short L_mis[MAX_K+1][2*MAX_K+1], uint8_t* quality, uint16_t* mis_n)
{
 
    _ASSERT(k < MAX_K);
    k = __min(MAX_K - 1, k); // enforce limit even in non-debug builds
    if (NULL == text) {
        // This happens when we're trying to read past the end of the genome.
        return -1;
    }
	
	uint8_t* quality_p = quality;
	uint16_t mis_tol = 0;
	_uint64 x = 0;
	_uint64 x_tmp = 0;
	_uint64 x_tmp_tmp = 0;
	
    const char* p = pattern;
    const char* t = text;
    int endl = __min(patternLen, textLen);
    const char* pend = pattern + endl; 
    while (p < pend) {
        
        x_tmp = *((_uint64*) p) ^ *((_uint64*) t);
		x = (x_tmp & (*((_uint64*)quality_p)));
		
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, endl);
			
			x_tmp_tmp = x_tmp & bv_64[zeroes];
			mis_tol += popcount_3((((x_tmp_tmp & low_bit_mask_lv) << 1) | (x_tmp_tmp & high_bit_mask_lv)));
			L_mis[0][MAX_K] = mis_tol;
			
            goto done1;
        }else	mis_tol += popcount_3((((x_tmp & low_bit_mask_lv) << 1) | (x_tmp & high_bit_mask_lv)));

        p += 8;
        t += 8;
		quality_p += 8;
    }
    L[0][MAX_K] = endl;
done1:
    if (L[0][MAX_K] == endl) {
        int result = (patternLen > endl ? patternLen - endl : 0); // Could need some deletions at the end
		
		(*mis_n) = mis_tol;
        return result;
    }
    int e, d;
    for ( e = 1; e <= k; e++) {
        // Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
        for ( d = 0; d != e+1; d = (d > 0 ? -d : -d+1)) {
            int best = L[e-1][MAX_K+d] + 1; // up
			mis_tol = L_mis[e-1][MAX_K+d];
			
            int left = L[e-1][MAX_K+d-1];
            if (left > best)
			{
				best = left;
				mis_tol = L_mis[e-1][MAX_K+d-1];
			}
                
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best)
			{
				best = right;
				mis_tol = L_mis[e-1][MAX_K+d+1];
			}
			
            const char* p = pattern + best;
            const char* t = (text + d) + best;
			quality_p = quality + best;
			
            if (*p == *t) { 
                int endl = __min(patternLen, textLen - d);
                const char* pend = pattern + endl;

                while (1) {
    
                    x_tmp = *((_uint64*) p) ^ *((_uint64*) t);
					x = x_tmp & (*((_uint64*)quality_p));
					
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, endl);
						
						x_tmp_tmp = x_tmp & bv_64[zeroes];
						mis_tol += popcount_3((((x_tmp_tmp & low_bit_mask_lv) << 1) | (x_tmp_tmp & high_bit_mask_lv)));
									
                        break;
                    }else	mis_tol += popcount_3((((x_tmp & low_bit_mask_lv) << 1) | (x_tmp & high_bit_mask_lv)));
					
                    p += 8;
					quality_p += 8;
					
                    if (p >= pend) {
                        best = endl;
                        break;
                    }
                    t += 8;
                }
            }

            if (best == patternLen) {
				(*mis_n) = mis_tol + e;
                return e;
            }
            L[e][MAX_K+d] = best;
			L_mis[e][MAX_K+d] = mis_tol;
        }
    }
    return -1;
}





int computeEditDistanceWithCigar_s_mis(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], uint8_t* quality)//bool useM, CigarFormat format, 
{
	_ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
	
	uint8_t* quality_p = quality;
	
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
		
		x &= (*((_uint64*)quality_p));
		
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
		quality_p += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
			
			quality_p = quality + best;
			
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
					
					x &= (*((_uint64*)quality_p));
					
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
					quality_p += 8;
					
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
            
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
				
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}
							accumulatedMs = 0;
						}
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	max_L = curD;

	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];

        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	accumulatedMs = L[0][MAX_K+0];

    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}
				accumulatedMs = 0;
			}
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - L[e - 1][MAX_K + max_L], 'S')) {//, format
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string

	
    //return -1;
	if(max_L < 0)	max_L = -max_L;
	
	return max_L;
}




int computeEditDistanceWithCigar_s_mis_left(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], uint8_t* quality, uint16_t* s_offset)//bool useM, CigarFormat format, 
{
    _ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
	
	uint8_t* quality_p = quality;
	
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
		
		x &= (*((_uint64*)quality_p));
		
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
		quality_p += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
			
			quality_p = quality + best;
			
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
					
					x &= (*((_uint64*)quality_p));
					
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
					quality_p += 8;
					
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}
							accumulatedMs = 0;
						}
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	max_L = curD;
	
	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];
        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	accumulatedMs = L[0][MAX_K+0];
	
    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}
				accumulatedMs = 0;
			}
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - L[e - 1][MAX_K + max_L], 'S')) {//, format
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string

	(*s_offset) = patternLen - L[e - 1][MAX_K + max_L];
    
	//return -1;
	if(max_L < 0)	max_L = -max_L;
	
	return max_L;
}

int computeEditDistanceWithCigar_s_nm_left(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], int* nm_score, uint16_t* s_offset)//bool useM, CigarFormat format, 
{
    _ASSERT(k < MAX_K);
	
	int m_num = 0;
	
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1];
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}else{
								m_num += accumulatedMs;
							}
							accumulatedMs = 0;
						}
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}else{
						m_num += accumulatedMs;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
				
				(*nm_score) = patternLen - m_num;
				
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	int curD_re = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	curD_re	= curD;
	//max_L = curD;
	
	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];
        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	accumulatedMs = L[0][MAX_K+0];
	
    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}else{
					m_num += accumulatedMs;
				}
				accumulatedMs = 0;
			}
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}else{
			m_num += accumulatedMs;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - max_L, 'S')) {//, format L[e - 1][MAX_K + max_L]
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string

	(*nm_score) = max_L - m_num;
	(*s_offset) = patternLen - max_L;
	
    //return -1;
	if(curD_re < 0)	curD_re = -curD_re;
	
	return curD_re;
}





int computeEditDistanceWithCigar_s_left(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], uint16_t* s_offset)//bool useM, CigarFormat format, 
{
    _ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}
							accumulatedMs = 0;
						}
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	max_L = curD;
	
	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];
        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	accumulatedMs = L[0][MAX_K+0];
	
    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}
				accumulatedMs = 0;
			}
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - L[e - 1][MAX_K + max_L], 'S')) {//, format
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	(*s_offset) = patternLen - L[e - 1][MAX_K + max_L];
	
    //return -1;
	if(max_L < 0)	max_L = -max_L;
	
	return max_L;
}



//for mapping quality

int computeEditDistanceWithCigar_s_mp(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], float* mp_subs, float* sub_t)//bool useM, CigarFormat format, 
{
    _ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
	uint16_t quality_offset = 0;
	uint16_t quality_i = 0;
	
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					
					int streakStart = 0;
					bool matching = (pattern[0] == text[0]);
					for (i = 0; i < end; i++) {
						bool newMatching = (pattern[i] == text[i]);
						if (newMatching != matching) {
							if(!matching)
							{
								//cal mapping quality
#ifdef	LV_MP
								for(quality_i = 0; quality_i < i - streakStart; quality_i++)
									(*sub_t) += mp_subs[quality_offset + quality_i];
#endif
							}
							quality_offset += i - streakStart;
							
							matching = newMatching;
							streakStart = i;
						}
					}
					
					// Write the last '=' or 'X' streak
					if (patternLen > streakStart) {
						if (!matching) {
							// Write out X's all the way to patternLen
								
							//cal mapping quality
#ifdef	LV_MP
							for(quality_i = 0; quality_i < patternLen - streakStart; quality_i++)		
								(*sub_t) += mp_subs[quality_offset + quality_i];	
#endif							
						} else {
							// Write out some ='s and then possibly X's if pattern is longer than text
							quality_offset += end - streakStart;
								
							if (patternLen > end) {
								
								//cal mapping quality
#ifdef	LV_MP
								for(quality_i = 0; quality_i < patternLen - end; quality_i++)		
									(*sub_t) += mp_subs[quality_offset + quality_i];	
#endif								
							}
						}
					}
						
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}

                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				if (L[0][MAX_K+0] > 0) {
					quality_offset += L[0][MAX_K+0];
						
				}
					
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					if((action == '=') || (action == 'I'))
						quality_offset += actionCount;
						
					if(action == 'X')
					{
						//cal mapping quality
#ifdef	LV_MP
						for(quality_i = 0; quality_i < actionCount; quality_i++)		
							(*sub_t) += mp_subs[quality_offset + quality_i];	
#endif						
					}
						
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}
							accumulatedMs = 0;
						}
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						
						quality_offset += backtraceMatched[curE];
						
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	max_L = curD;
	
	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];
        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	if (L[0][MAX_K+0] > 0) {
		quality_offset += L[0][MAX_K+0];
	}
	
	accumulatedMs = L[0][MAX_K+0];
	
    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		
		if((action == '=') || (action == 'I'))
			quality_offset += actionCount;
						
		if(action == 'X')
		{
			//cal mapping quality
#ifdef	LV_MP
			for(quality_i = 0; quality_i < actionCount; quality_i++)		
				(*sub_t) += mp_subs[quality_offset + quality_i];
#endif	
			quality_offset += actionCount;
		}
	
				
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}
				accumulatedMs = 0;
			}
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			
			quality_offset += backtraceMatched[curE];
			
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - L[e - 1][MAX_K + max_L], 'S')) {//, format
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string

	
    //return -1;
	if(max_L < 0)	max_L = -max_L;
	
	return	max_L;
}



int computeEditDistanceWithCigar_s_mis_mp(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], uint8_t* quality, float* mp_subs, float* sub_t)//bool useM, CigarFormat format, 
{
    _ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
	
	uint8_t* quality_p = quality;
	
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
		
		x &= (*((_uint64*)quality_p));
		
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
		quality_p += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
	uint16_t quality_offset = 0;
	uint16_t quality_i = 0;
	 
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
			
			quality_p = quality + best;
			
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
					
					x &= (*((_uint64*)quality_p));
					
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
					quality_p += 8;
					
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					int streakStart = 0;
					bool matching = (pattern[0] == text[0]);
					for (i = 0; i < end; i++) {
						bool newMatching = (pattern[i] == text[i]);
						if (newMatching != matching) {
							
							if(!matching)
							{
								//cal mapping quality
#ifdef	LV_MP
								for(quality_i = 0; quality_i < i - streakStart; quality_i++)
									(*sub_t) += mp_subs[quality_offset + quality_i];
#endif
							}
							quality_offset += i - streakStart;

							matching = newMatching;
							streakStart = i;
						}
					}
					
					// Write the last '=' or 'X' streak
					if (patternLen > streakStart) {
						if (!matching) {
							// Write out X's all the way to patternLen
							//cal mapping quality
#ifdef	LV_MP
							for(quality_i = 0; quality_i < patternLen - streakStart; quality_i++)		
								(*sub_t) += mp_subs[quality_offset + quality_i];			
#endif
						} else {
							// Write out some ='s and then possibly X's if pattern is longer than text
							quality_offset += end - streakStart;
							
							if (patternLen > end) {
								//cal mapping quality
#ifdef	LV_MP
								for(quality_i = 0; quality_i < patternLen - end; quality_i++)		
									(*sub_t) += mp_subs[quality_offset + quality_i];	
#endif								
							}
						}
					}

					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				if (L[0][MAX_K+0] > 0) {
					quality_offset += L[0][MAX_K+0];
						
				}
				
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					if((action == '=') || (action == 'I'))
						quality_offset += actionCount;
						
					if(action == 'X')
					{
						//cal mapping quality
#ifdef	LV_MP
						for(quality_i = 0; quality_i < actionCount; quality_i++)		
							(*sub_t) += mp_subs[quality_offset + quality_i];	
#endif		
					}
					
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}
							accumulatedMs = 0;
						}
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						
						quality_offset += backtraceMatched[curE];
						
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	max_L = curD;
	
	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];
        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	if (L[0][MAX_K+0] > 0) {
		quality_offset += L[0][MAX_K+0];
	}
	
	accumulatedMs = L[0][MAX_K+0];
	
    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		if((action == '=') || (action == 'I'))
			quality_offset += actionCount;
						
		if(action == 'X')
		{
			//cal mapping quality
#ifdef	LV_MP
			for(quality_i = 0; quality_i < actionCount; quality_i++)		
				(*sub_t) += mp_subs[quality_offset + quality_i];
#endif	
			quality_offset += actionCount;
		}
		
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}
				accumulatedMs = 0;
			}
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			quality_offset += backtraceMatched[curE];
			
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - L[e - 1][MAX_K + max_L], 'S')) {//, format
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string

	
    //return -1;
	if(max_L < 0)	max_L = -max_L;
	
	return max_L;
}




int computeEditDistanceWithCigar_s_mis_left_mp(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, 
    short L[MAX_K+1][2*MAX_K+1], uint8_t* quality, uint16_t* s_offset, float* mp_subs, float* sub_t)//bool useM, CigarFormat format, 
{
    _ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.
#ifdef LV_INI
    char A[MAX_K+1][2*MAX_K+1] = {0};
    char backtraceAction[MAX_K+1] = {0};
    int backtraceMatched[MAX_K+1] = {0};
    int backtraceD[MAX_K+1] = {0};
#else
	char A[MAX_K+1][2*MAX_K+1];
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
#endif
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
	
	uint8_t* quality_p = quality;
	
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
		
		x &= (*((_uint64*)quality_p));
		
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
		quality_p += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		/*
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=')) {//, format
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
					return -2;
				}
			}
		}
		*/
		
		if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
				return -2;
			}
		
        return 0;
    }
    int e;
	int best;
	int best_i;
	int curE;
	int curD;
	
	uint16_t quality_offset = 0;
	uint16_t quality_i = 0;
	
    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        int d;
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
			
			quality_p = quality + best;
			
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
					
					x &= (*((_uint64*)quality_p));
					
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
					quality_p += 8;
					
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.
                
                int straightMismatches = 0;
                int i;
                for (i = 0; i < end; i++) {
                    if (pattern[i] != text[i]) {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e) {
					
					/*
                    // We can match with no indels; let's do that
					if (useM) {
						//
						// No inserts or deletes, and with useM equal and SNP look the same, so just
						// emit a simple string.
						//
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'))) {//, format
									return -2;
								}
								matching = newMatching;
								streakStart = i;
							}
						}
					
						// Write the last '=' or 'X' streak
						if (patternLen > streakStart) {
							if (!matching) {
								// Write out X's all the way to patternLen
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X')) {//, format
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=')) {//, format
									return -2;
								}
								if (patternLen > end) {
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X')) {//, format
										return -2;
									}
								}
							}
						}
					}
					*/
					int streakStart = 0;
					bool matching = (pattern[0] == text[0]);
					for (i = 0; i < end; i++) {
						bool newMatching = (pattern[i] == text[i]);
						if (newMatching != matching) {
							if(!matching)
							{
								//cal mapping quality
#ifdef	LV_MP
								for(quality_i = 0; quality_i < i - streakStart; quality_i++)
									(*sub_t) += mp_subs[quality_offset + quality_i];
#endif
							}
							quality_offset += i - streakStart;
							
							matching = newMatching;
							streakStart = i;
						}
					}
					
					// Write the last '=' or 'X' streak
					if (patternLen > streakStart) {
						if (!matching) {
							// Write out X's all the way to patternLen
#ifdef	LV_MP	
							//cal mapping quality
							for(quality_i = 0; quality_i < patternLen - streakStart; quality_i++)		
								(*sub_t) += mp_subs[quality_offset + quality_i];
#endif							
						} else {
							// Write out some ='s and then possibly X's if pattern is longer than text
							quality_offset += end - streakStart;
								
							if (patternLen > end) {
#ifdef	LV_MP
								//cal mapping quality
								for(quality_i = 0; quality_i < patternLen - end; quality_i++)		
									(*sub_t) += mp_subs[quality_offset + quality_i];
#endif								
							}
						}
					}
					
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M')) {//, format
							return -2;
						}
					
                    return e;
                }
                
#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                int ee;
                for (ee = 0; ee <= e; ee++) {
                    int dd;
                    for (dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                curD = d;
                //int curE;
                for (curE = e; curE >= 1; curE--) {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I') {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    } else if (backtraceAction[curE] == 'D') {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    } else { // backtraceAction[curE] == 'X'
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                        backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				/*
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
							return -2;
						}
					}
				}
				*/
				if (L[0][MAX_K+0] > 0) {
					quality_offset += L[0][MAX_K+0];
						
				}
				
				accumulatedMs = L[0][MAX_K+0];

                curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					
					/*
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
									return -2;
								}
								accumulatedMs = 0;
							}
							if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
								return -2;
							}
						}
					} else {
						if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					*/
					
					if((action == '=') || (action == 'I'))
						quality_offset += actionCount;
						
					if(action == 'X')
					{
#ifdef	LV_MP
						//cal mapping quality
						for(quality_i = 0; quality_i < actionCount; quality_i++)		
							(*sub_t) += mp_subs[quality_offset + quality_i];	
#endif				
					}
					
					if (action == '=' || action == 'X') {
						accumulatedMs += actionCount;
					} else {
						if (accumulatedMs != 0) {
							if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
								return -2;
							}
							accumulatedMs = 0;
						}
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
							return -2;
						}
					}
					
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						/*
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
							if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
								return -2;
							}
						}
						*/
						quality_offset += backtraceMatched[curE];
						
						accumulatedMs += backtraceMatched[curE];
                    }
                    curE++;
                }
				if (1 && accumulatedMs != 0) {
					//
					// Write out the trailing Ms.
					//
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
				}
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    //*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	
	int max_L = 0;
	for(best_i = 1 - e; best_i < e; best_i++)
		if(L[e - 1][MAX_K + best_i] > max_L)
		{
			max_L = L[e - 1][MAX_K + best_i];
			curD = best_i;
		}
		
	max_L = curD;
	
	//for cigar

    //int curE;
    for (curE = e - 1; curE >= 1; curE--) {
        backtraceAction[curE] = A[curE][MAX_K+curD];
        if (backtraceAction[curE] == 'I') {
            backtraceD[curE] = curD + 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
        } else if (backtraceAction[curE] == 'D') {
            backtraceD[curE] = curD - 1;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
        } else { // backtraceAction[curE] == 'X'
            backtraceD[curE] = curD;
            backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
        }
        curD = backtraceD[curE];
    }

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	/*
	if (useM) {
		accumulatedMs = L[0][MAX_K+0];
	} else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K+0] > 0) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=')) {//, format
				return -2;
			}
		}
	}
	*/
	if (L[0][MAX_K+0] > 0) {
		quality_offset += L[0][MAX_K+0];
	}
	
	accumulatedMs = L[0][MAX_K+0];
	
    curE = 1;
    while (curE <= e - 1) {
        // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
        char action = backtraceAction[curE];
        int actionCount = 1;
        while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
            actionCount++;
            curE++;
        }
		
		/*
		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			} else {
				if (accumulatedMs != 0) {
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
						return -2;
					}
					accumulatedMs = 0;
				}
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
						return -2;
				}
			}
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
				return -2;
			}
		}
		*/
		if((action == '=') || (action == 'I'))
			quality_offset += actionCount;
						
		if(action == 'X')
		{
			//cal mapping quality
#ifdef	LV_MP
			for(quality_i = 0; quality_i < actionCount; quality_i++)		
				(*sub_t) += mp_subs[quality_offset + quality_i];
#endif			
			quality_offset += actionCount;
		}
		
		if (action == '=' || action == 'X') {
			accumulatedMs += actionCount;
		} else {
			if (accumulatedMs != 0) {
				if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
					return -2;
				}
				accumulatedMs = 0;
			}
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action)) {//, format
					return -2;
			}
		}
			
        // Next, write out ='s for the exact match
        if (backtraceMatched[curE] > 0) {
			/*
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			} else {
				if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=')) {//, format
					return -2;
				}
			}
			*/
			
			quality_offset += backtraceMatched[curE];
			
			accumulatedMs += backtraceMatched[curE];
        }
        curE++;
	}
	if (1 && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M')) {//, format
			return -2;
		}
	}
	
	if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - L[e - 1][MAX_K + max_L], 'S')) {//, format
		return -2;
	}
		
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string

	(*s_offset) = patternLen - L[e - 1][MAX_K + max_L];
    
	//return -1;
	if(max_L < 0)	max_L = -max_L;
	
	return max_L;
}