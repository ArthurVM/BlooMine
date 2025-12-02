#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

// constants for affine alignment scoring
const double hit = 10.0;          // the value to increment the score by for a hit event
const double gap_open = 15.0;     // the value to penalise the score by for a gap open event
const double neg = 7.0;           // the value to penalise the score by for a gap extension event

#endif /* CONSTANTS_HPP */