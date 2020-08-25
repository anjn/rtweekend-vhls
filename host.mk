root    ?= $(shell readlink -f $(dir $(firstword $(MAKEFILE_LIST))))
cppmake ?= $(shell readlink -f $(dir $(lastword  $(MAKEFILE_LIST))))

build_dir ?= .
bin_dir   ?= $(build_dir)/bin
lib_dir   ?= $(build_dir)/lib
obj_dir   ?= $(build_dir)/obj

cpp_ext   ?= .cpp

target_object_files = $(patsubst %, $(obj_dir)/%.o, $(targets))
object_files = $(patsubst %, $(obj_dir)/%.o, $(objects))
all_object_files = $(target_object_files) $(object_files)
depend_files = $(patsubst %.o,%.d,$(all_object_files))
bin_files = $(addprefix $(bin_dir)/, $(targets))
lib_files = $(addprefix $(lib_dir)/, $(libraries))

.DEFAULT_GOAL := all
all: $(depend_targets) submodules targets

submodules: $(submodules)

targets: $(bin_files) $(lib_files)

# include projects
link_libs := $(link_libs)
include $(addsuffix /include.mk, $(include_prjs))

.PHONY: all submodules targets clean $(clean_targets)
.SECONDARY: $(all_object_files)

RANLIB    ?= ranlib
CXXFLAGS  += -Wall -Wextra
CXXFLAGS  += -Wno-unused-parameter -Wno-sign-compare -Wno-missing-field-initializers
CXXFLAGS  += -Wformat=2 -Wstrict-aliasing=2 -Wdisabled-optimization
CXXFLAGS  += -Wfloat-equal -Wpointer-arith -Wcast-align -Wredundant-decls
CXXFLAGS  += $(addprefix -I, $(include_dirs))
LDFLAGS   += 
MAKEFLAGS += --no-print-directory

ifdef DEBUG
	CXXFLAGS += -O0 -g -gdwarf-3 -fno-inline
else
	CXXFLAGS += -O3
endif

$(bin_dir)/%: $(obj_dir)/%.o $(object_files) $(link_libs) $(bin_depends)
	@mkdir -p $(bin_dir)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(lib_dir)/lib%.so: $(object_files) $(lib_depends)
	@mkdir -p $(lib_dir)
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

$(lib_dir)/lib%.a: $(object_files) $(lib_depends)
	@mkdir -p $(lib_dir)
	$(AR) cr $@ $^
	$(RANLIB) $@

$(obj_dir)/%.o: %$(cpp_ext) $(obj_depends)
	@mkdir -p $(dir $@)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean: $(clean_targets)
	rm -f $(bin_files) $(lib_files) $(all_object_files) $(depend_files)
	test $(bin_dir) != "." && rm -rf $(bin_dir) || true
	test $(lib_dir) != "." && rm -rf $(lib_dir) || true
	test $(obj_dir) != "." && rm -rf $(obj_dir) || true

# depend
$(obj_dir)/%.d: %$(cpp_ext)
	@mkdir -p $(dir $@)
	@rm -f $@
	@$(CXX) -MM $(CXXFLAGS) $< | \
		sed -e 's/$(subst /,\/,$*)\.o[ :]*/$(subst /,\/,$(dir $@))$(subst /,\/,$*).o $(subst /,\/,$(dir $@))$(subst /,\/,$*).d : $$\(wildcard /' | tr -d '\\\n' | tr -d '\\\r\n' > $@ || true; \
		[ -s $@ ] && echo ')' >> $@; \
		[ -s $@ ] || rm -f $@

-include $(depend_files)
