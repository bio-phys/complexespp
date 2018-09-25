
#include <altivec.h>

int main() {
    {
        __vector double res0d;
        __vector double res1d;
        __vector double res2d = vec_add(res0d, res1d);

        res2d = vec_abs(res0d);
        res2d = vec_rsqrt(res0d);
    }
    {
        __vector float res0;
        __vector float res1;
        __vector float res2 = vec_add(res0, res1);

        res2 = vec_abs(res0);
        res2 = vec_rsqrt(res0);
    }
    return 0;
}
