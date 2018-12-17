#include "nrg.h"

extern float good,cons;
extern float del2,delta;
extern int idel;

void Ajustadelta(void)
{
    float fail;

    fail=1.0F - good/(cons * 4);
    if (fail > 0.6F  &&  delta > 0.01F )
	delta *= 0.9F;
    else if (fail < 0.6F  && delta < 4.0F )
	delta *= 1.1F;
    idel=(int) (delta/twopioverNest);
    del2= 2.0F *( (float) idel) + 1.0F;
}
