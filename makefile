fc      = ifort
fcflags = -g -O3 -mavx
src     = gauss_quad1d.F90
bin     = a.out

$(bin): $(src)
	$(fc) $(fcflags) $^ -o $@

all: $(bin)

clean:
	rm -f $(bin) *.mod
