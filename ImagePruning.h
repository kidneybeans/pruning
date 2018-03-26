#ifndef _IMAGEPRUNING_H
#define _IMAGEPRUNING_H

#include <string>
#include <vector>
#include <stdio.h>
#include <fstream>
#include "ip.h"
#include "graphics.h"
#include "public.h"
#include "drwnBase.h"
#include "zlib.h"
#include "png.h"
#include <stdlib.h>

#include "colorseg.h"

using namespace std;

typedef struct  tag_WebImage{

	string imagepathes;
	int idx;
	int label;

}WebImage;

typedef struct tag_PruImage{

	CImage *binary_img;
	int line_num;
	rect *textlines;

#ifdef DEBUGERROR
	int etype;
#endif

public: 
	tag_PruImage():binary_img(NULL), line_num(0), textlines(NULL)
#ifdef DEBUGERROR
		,etype(0)
#endif
	{}

}PruImage;

void buildDatabase(string datasetdir, vector<vector<WebImage>> &dir_array);
void parseImagePruning(vector<vector<WebImage>> &dir_array, vector<vector<WebImage>> &pru_array);
void imageProcess(CImage &image, PruImage &pruimg, int &label);
void clearOutImage(PruImage &pruimg);
void evaluation(vector<vector<WebImage>> dir_array, vector<vector<WebImage>> pru_array);
void imageResize(CImage &image, CImage &dsimage);
void resizeImageBilinear(const uchar* src, const int src_width, const int src_height, const int src_step, uchar* dst, const int dst_width, 
	const int dst_height, const int dst_step, const int channels);
int convert(CImage &cImg, char* inname);
inline int findStr(char *str, char *substr)
{
	int  n;
	char  *p, *r;
	n = 0;
	while (*str)
	{
		p = str;
		r = substr;
		while (*r){
			if (*r == *p)
			{
				r++;
				p++;
			}
			else
			{
				break;
			}
		}
		if (*r == '\0')
			n++;
		str++;
	}
	return n;
}


static inline int round(double value)
{
	int t;
	__asm
	{
		fld value;
		fistp t;
	}
	return t;
}
static inline short saturate_cast_short(int v)
{
	return (short)((unsigned)(v - SHRT_MIN) <= (unsigned)USHRT_MAX ? v : v > 0 ? SHRT_MAX : SHRT_MIN);
}
static inline short saturate_cast_short(float v)
{ 
	int iv = round(v); 
	return saturate_cast_short(iv); 
}
// helper function
static inline int clip(int x, int a, int b)
{
	return x >= a ? (x < b ? x : b-1) : a;
}

#endif