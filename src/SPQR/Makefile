# compile just the C-callable library
# need to make a target whose name is not lib because on Windows "Lib" == "lib"

.PHONY: library clean

library:
	( cd Lib ; $(MAKE) lib )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Lib ; $(MAKE) clean )
