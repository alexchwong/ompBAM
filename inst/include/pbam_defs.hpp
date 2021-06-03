#ifndef _pbam_defs
#define _pbam_defs

// Definitions

static const int bamGzipHeadLength = 16;  // +2 a uint16 with the full block length.
static const char bamGzipHead[bamGzipHeadLength+1] = 
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00";
static const int bamEOFlength = 28;
static const char bamEOF[bamEOFlength+1] =
		"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";

static const int magiclength = 4;
static const char magicstring[magiclength+1] = "\x42\x41\x4d\x01";

#endif