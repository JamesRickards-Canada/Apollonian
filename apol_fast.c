/*Methods to quickly compute Apollonian stuff, and supporting things.*/

/*INCLUSIONS*/
#include <pari/pari.h>
#include "apol.h"




/*SECTION 1: MISSING CURVATURES*/


/*MOVE THIS TO APOL.C*/

/*Returns the "type" of v, i.e. the number of resiudes modulo 24 and the smallest residue coprime to 6, which uniquely identifies it.*/
GEN
apol_type(GEN v)
{
  pari_sp av = avma;
  GEN m24 = apol_mod24(v);
  long second = m24[2];/*Uniquely identified by the second element*/
  set_avma(av);
  switch (second) {
    case 1: return mkvec2s(6, 1);
	case 3: return mkvec2s(8, 11);
	case 4: return mkvec2s(6, 13);
	case 5: return mkvec2s(6, 5);
	case 6: return mkvec2s(8, 7);
	case 8: return mkvec2s(6, 17);
  }
  pari_err(e_MISC, "We didn't find one of the 6 possible admissible sets, are you sure you inputted a primitive Apollonian circle packing?");
  return gen_0;
}



