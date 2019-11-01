
target=$1

# ======= DEFS / SETUP

BUILDDIR="./Build" 

IFLAGS="-I./math -I./GUI"
LFLAGS="-L/usr/lib/x86_64-linux-gnu  -lGL -lSDL2 -lm"
#          /usr/lib/x86_64-linux-gnu/

#CFLAGS=" -std=c99 -Og"
#CFLAGS="-std=c99 -Og -Wall"
CFLAGS="-std=c99 -Og -Wall -Wno-missing-braces -Werror=return-type"
#CFLAGS="-E"
#CFLAGS=" -std=c99 -O3 -ffast-math -ftree-vectorize" 
#CFLAGS=" -std=c99 -O2 -march=native -mtune=native"
#CFLAGS=" -std=c99 -O3 -march=native -mtune=native"
#CFLAGS="-std=c99 -Ofast -march=native -mtune=native"

# ====== Aux Defs

BINAME=$BUILDDIR/$target.x

# ====== Main

mkdir $BUILDDIR

rm $BINAME

#gcc $IFLAGS -E $target.c
gcc $CFLAGS $IFLAGS $target.c -o $BINAME $LFLAGS
$BINAME