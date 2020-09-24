#!/usr/bin/env python3
import os

def fix_file(filename, to_replace):
    lines = []
    with open(filename) as f:
        for line in f:
            if line in to_replace:
                print(f"Replacing in {filename}:\n", line, to_replace[line], sep="")
                line = to_replace[line]
            lines.append(line)

    with open(filename, "w") as f:
        print(*lines, sep="", file=f)


filename = os.path.join("mccortex", "src", "global", "util.h")
to_replace = {
    "const uint8_t rev_nibble_table[16];\n": "extern const uint8_t rev_nibble_table[16];\n"
}
fix_file(filename, to_replace)


filename = os.path.join("mccortex", "libs", "xxHash", "xxhsum.c")
to_replace = {
    "#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(_WIN32) || defined(__CYGWIN__)\n": "#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(_WIN32)\n"
}
fix_file(filename, to_replace)


filename = os.path.join("mccortex", "libs", "bcftools", "Makefile")
to_replace = {
   "DYNAMIC_FLAGS = -rdynamic\n": "DYNAMIC_FLAGS = \n",
   "ALL_LIBS     = -lz -ldl $(LIBS)\n": "ALL_LIBS     = -lz $(LIBS)\n",
   "PLUGINS_ENABLED = yes\n": "PLUGINS_ENABLED = no\n"
}
fix_file(filename, to_replace)


filename = os.path.join("mccortex", "libs", "seq_file", "Makefile")
to_replace = {
   "LINKING=$(HTSARGS) -lpthread -lz\n": "LINKING=$(HTSARGS) -lpthread -lws2_32 -lz -llzma -lcurl -lbz2\n",
}
fix_file(filename, to_replace)


filename = os.path.join("mccortex", "Makefile")
to_replace = {
   "LINK=-lpthread -lz -lm\n": "LINK=-lpthread -lz -lm -lws2_32 -lz -llzma -lcurl -lbz2\n",
}
fix_file(filename, to_replace)


filename = os.path.join("mccortex", "libs", "seq_file", "benchmarks", "Makefile")
to_replace = {
   "LINKING=$(HTSARGS) $(SAMLINK) -lpthread -lz\n": "LINKING=$(HTSARGS) $(SAMLINK) -lpthread -lz -lws2_32 -lz -llzma -lcurl -lbz2\n",
}
fix_file(filename, to_replace)


filename = os.path.join("mccortex", "libs", "seq_file", "dev", "Makefile")
to_replace = {
       "LINKING=$(HTSARGS) $(SAMLINK) -lpthread -lz\n": "LINKING=$(HTSARGS) $(SAMLINK) -lpthread -lz -lws2_32 -lz -llzma -lcurl -lbz2\n",
}
fix_file(filename, to_replace)


filename = os.path.join("mccortex", "libs", "readsim", "Makefile")
to_replace = {
       "\tLIBS=-lpthread -lz -lm\n": "\tLIBS=-lpthread -lz -lm -lws2_32 -lz -llzma -lcurl -lbz2\n",
}
fix_file(filename, to_replace)


filename = os.path.join("mccortex", "libs", "vcf-slim", "Makefile")
to_replace = {
       "LIBS=-lz -lm -lpthread\n": "LIBS=-lpthread -lz -lm -lws2_32 -lz -llzma -lcurl -lbz2\n",
}
fix_file(filename, to_replace)
