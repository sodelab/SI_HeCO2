#ifndef POLY_2B_H
#define POLY_2B_H

//
// reduced (no bonds breaking) permutational symmetry
// including only 2B terms
//
//    x[0] = @VAR@-|x-intra-CO|(Cb, Ob1);
//    x[1] = @VAR@-|x-intra-CO|(Cb, Ob2);
//    x[2] = @VAR@-|x-intra-OO|(Ob1, Ob2);
//    x[3] = @VAR@-|x-ArC|(Ara, Cb);
//    x[4] = @VAR@-|x-ArO|(Ara, Ob1);
//    x[5] = @VAR@-|x-ArO|(Ara, Ob2);
//
namespace x2b {

void eval(const double x[6], double*);

}

#endif // POLY_2B_H
