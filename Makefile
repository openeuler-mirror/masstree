
ifeq ($(strip $(MEMMGR)), )
  MEMMGR = -ljemalloc
endif
ifneq ($(strip $(KEYSWAP)), )
  CPPFLAGS += -DKEYSWAP
endif
ifneq ($(strip $(NOPREFETCH)), )
  CPPFLAGS += -DNOPREFETCH
endif
ifneq ($(strip $(NOSUPERPAGE)), )
  CPPFLAGS += -DNOSUPERPAGE
endif
LIBS = -lnuma -lpthread -lm
LDFLAGS = -fstack-protector-all -z,now

CXXFLAGS = -g -W -O3 -std=gnu++11 -Wextra -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC -g -Werror -DNDEBUG -D_FORTIFY_SOURCE=2 -Wall -Werror -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unused-function -Wno-strict-aliasing  -faligned-new -Wwrite-strings -Wcast-align -Wreturn-type -Wpointer-arith -Wlogical-op -Waddress -Wsizeof-pointer-memaccess -Winit-self -fno-exceptions -fno-rtti -Wnon-virtual-dtor -Wno-missing-field-initializers

ifeq ($(shell uname -p),aarch64)
	CXXFLAGS += -march=armv8-a+crc
else
	CXXFLAGS += -mcx16
endif

all: libmasstree.so

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(MEMMGR) $(LDFLAGS) $(LIBS) -c -o $@ $<

libmasstree.so: straccum.o string.o kvthread.o
	@rm -f $@
	$(CXX) -fstack-protector-all -Wl,-z,relro,-z,now -shared $^ -o $@

clean:
	rm -f *.o libmasstree.so

.PHONY: clean all
