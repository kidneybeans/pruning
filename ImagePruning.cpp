#include "stdafx.h"
#include "time.h"
#include "ImagePruning.h"

//build database: record the ground truth in 'dir_array'
void buildDatabase(string datasetdir, vector<vector<WebImage>> &dir_array)
{
	int idx = 0; 

	set<const char*> imagefiletypes;
	imagefiletypes.insert("bmp");
	imagefiletypes.insert("jpeg");
	imagefiletypes.insert("jpg");
	imagefiletypes.insert("tiff");
	imagefiletypes.insert("tif");
	imagefiletypes.insert("png");
	imagefiletypes.insert("gif");

	vector<WebImage> image_array, type_array;
	vector<string> imagesinaclass, typeinclass;
	vector<string> traningdirs=drwnDirectoryListing(datasetdir.c_str());
	dir_array.resize(traningdirs.size());

	for(unsigned int i=0;i<traningdirs.size();i++){
		imagesinaclass = drwnDirectoryListing(traningdirs[i].c_str(),imagefiletypes);

		idx=0;
		image_array.resize(imagesinaclass.size());
		for(unsigned int j=0;j<imagesinaclass.size();j++){
			if(findStr((char*)strWithoutExt(strFilename(imagesinaclass[j])).c_str(),"address")>0)
				image_array[j].label=0;
			else image_array[j].label=2;
			image_array[j].imagepathes=imagesinaclass[j];
			image_array[j].idx = idx++;
		}

		dir_array[i]=image_array;
	}

	for(unsigned int i=0;i<traningdirs.size();i++){
		string temp = traningdirs[i] + "\\type1";
		typeinclass = drwnDirectoryListing(temp.c_str());

		type_array.resize(typeinclass.size());
		for(unsigned int j=0;j<typeinclass.size();j++){
			for(unsigned int k=0;k<dir_array[i].size();k++){
				if(strcmp(strFilename(typeinclass[j]).c_str(), strFilename(dir_array[i][k].imagepathes).c_str())==0){
					if(dir_array[i][k].label==2) dir_array[i][k].label=1;
					break;
				}
			}
		}
	}
}

//get the pruning result: record the result in 'pru_array'
void parseImagePruning(vector<vector<WebImage>> &dir_array, vector<vector<WebImage>> &pru_array)
{
	//Parse Image
	int minimagesize=20;
	int maximagesize=200;
	int maxsize = 1000;//not support
	int label = 0;
	bool bresize = false;

	double times = 0;
	int numsamples = 0;
	for(unsigned int i=0; i<dir_array.size();i++){
		char buffer[20] = "";
		_itoa_s(i+1, buffer, 10);
		string s(buffer);

		for(unsigned int j=0;j<dir_array[i].size();j++){
			CImage image, dsimage;
			PruImage pruimg;

			bresize = false;

#ifdef DEBUG
			//if(strcmp("address.gif", strFilename(dir_array[i][j].imagepathes).c_str())!=0) //content_con1img4
			if(findStr((char*)strWithoutExt(strFilename(dir_array[i][j].imagepathes)).c_str(),"359b033b5bb5c9ea396d1d71d539b6003bf3b38"/*"footer_co"*/)<=0)
				continue;		
#endif
			//if the image is a 'png' file, then use this function to get the image data
			if(strcmp((char*)strExtension(strFilename(dir_array[i][j].imagepathes)).c_str(),"png")==0){
				convert(image, (char*)dir_array[i][j].imagepathes.c_str());
				//image.WriteImage("test.bmp");
			}
			else
				image.ReadImage((char*)dir_array[i][j].imagepathes.c_str(), true);

			//too small images and too large (height) images are not considered as address images
			if(image.GetWidth()<minimagesize&&image.GetHeight()<minimagesize){
				pru_array[i][j].label = 2;
				//very important: we do not support to large images
				dir_array[i][j].label = 2;
				goto SAVE;
			}
			if(/*image.GetWidth()>=maxsize||*/image.GetHeight()>=maxsize){
				pru_array[i][j].label = 2;
				//very important: we do not support to large images
				dir_array[i][j].label = 2;
				goto SAVE;
			}

			if(image.GetWidth()>maximagesize && image.GetHeight()>maximagesize){
				float ratio=((image.GetWidth()<image.GetHeight())? image.GetWidth():image.GetHeight())/(float(maximagesize));
				int w1 = (int)(image.GetWidth()/ratio), h1 = (int)(image.GetHeight()/ratio);
				dsimage.CreateImage(w1, h1, true);
				imageResize(image, dsimage);
				bresize = true;
			}

			//Process
			clock_t start,finish;
			double totaltime;
			start=clock();

			//get the image label and output text line informations
			if (bresize) imageProcess(dsimage, pruimg, label);
			else imageProcess(image, pruimg, label);

#ifdef DEBUGERROR
			int etype = 0;
			etype = pruimg.etype;
#endif
			//do something to pruimg then delete
			if(pruimg.binary_img) clearOutImage(pruimg);

			finish=clock();
			totaltime=(double)(finish-start);

			times += totaltime;
			numsamples ++;
			cout<<"\n processing time of "<< strFilename(dir_array[i][j].imagepathes).c_str() << " is "<<totaltime<<" milliseconds"<<endl;

			//Get the output array
			pru_array[i][j].label = label;

SAVE:
			//for debug
			char buffer1[20] = "", buffer2[20] = "";
			_itoa_s(dir_array[i][j].label, buffer1, 10);
			_itoa_s(pru_array[i][j].label, buffer2, 10); 
			string s1(buffer1), s2(buffer2); 

#ifdef DEBUGERROR
			char buffer3[20] = "";
			_itoa_s(etype, buffer3, 10);
			string s3(buffer3);
			string tmp = "..\\output\\" + s + "____" + strWithoutExt(strFilename(dir_array[i][j].imagepathes)) + "____" + s1 + "_" + s2 + "_" + s3 + ".bmp";
#else
			//string tmp = "..\\output\\" + strWithoutExt(strFilename(dir_array[i][j].imagepathes)) + "____" + s1 + "_" + s2 + "." + strExtension(dir_array[i][j].imagepathes);
			string tmp = "..\\output\\" + s + "____" + strWithoutExt(strFilename(dir_array[i][j].imagepathes)) + "____" + s1 + "_" + s2 + ".bmp";
#endif

			if(bresize) dsimage.WriteImage((char*)tmp.c_str());
			else image.WriteImage((char*)tmp.c_str());

			if(image.m_Data) {
				image.DeleteImage();
				image.m_Data = 0;
			}
			if(dsimage.m_Data){
				dsimage.DeleteImage();
				dsimage.m_Data = 0;
			}
			if(pruimg.binary_img) clearOutImage(pruimg);
		}
	}

	printf("\n average time on samples is: %.1f\n", times/numsamples);
}

void clearOutImage(PruImage &pruimg)
{
	if(pruimg.binary_img){
		if(pruimg.binary_img->m_Data){
			pruimg.binary_img->DeleteImage(); pruimg.binary_img->m_Data =0;
		}
		delete pruimg.binary_img; pruimg.binary_img = NULL;
	}

	if(pruimg.textlines){
		free(pruimg.textlines); pruimg.textlines = NULL;
	}
}

//get the image label and output text line informations
void imageProcess(CImage &image, PruImage &pruimg, int &label)
{
	//get the image label
	colorSeg Seg(image, false);
	label = Seg.getLabel();

	//get the output text line informations
	pruimg.line_num = Seg.getLinesNum();
	pruimg.textlines = Seg.getLines();
	pruimg.binary_img = Seg.getLblImg();

#ifdef DEBUGERROR
	pruimg.etype = Seg.getErrorType();
#endif
}

//evaluate the pruning result by compare the 'dir_array' and 'pru_array'
void evaluation(vector<vector<WebImage>> dir_array, vector<vector<WebImage>> pru_array)
{
	float recall = 0.f, precision = 0.f, totalrecall = 0.f;
	int t0 = 0, t1 = 0, t2 = 0;
	int det0 = 0, det1 = 0, det2 = 0;
	int tdet0 = 0, tdet1 = 0, tdet2 = 0;

	for(unsigned int i=0; i<dir_array.size();i++){
		for(unsigned int j=0;j<dir_array[i].size();j++){
			int tmp = dir_array[i][j].label;
			if(tmp == 0) t0++;
			else if(tmp ==1) t1++;
			else t2++;
		}
	}

	for(unsigned int i=0; i<dir_array.size();i++){
		for(unsigned int j=0;j<dir_array[i].size();j++){
			int truelabel = dir_array[i][j].label;
			int tmp = pru_array[i][j].label;
			if(tmp==0) { 
				if(truelabel==0) tdet0++;
			}
			if(tmp==1) {
				if(truelabel==0) tdet1++;
				if(truelabel==1) tdet1++;
				det1++;
			}
			if(tmp==2) {
				if(truelabel==2) tdet2++;
				det2++;
			}
		}
	}

	totalrecall = (float)tdet0/t0*100;
	recall = (float)tdet1/(t0+t1)*100;
	precision = (float)tdet1/det1*100;
	printf("Recall = %.1f%, Recall = %.1f%, Precision = %.4f%\n", totalrecall, recall, precision);
 	printf("Recall: det1 = %d, gt = %d, Precision: det1 = %d, gt = %d\n", tdet1, t0+t1, tdet1, det1);
}

void imageResize(CImage &image, CImage &dsimage)
{
	if(0){
		int ii, jj, ind, ind2, i;

		int w = image.GetWidth(), h = image.GetHeight();
		int w1 = dsimage.GetWidth(), h1 = dsimage.GetHeight();

		float f = 1.0f * w / w1;
		int *pbuf = new int[max(w1, h1)];

		for (ii = 0; ii < max(w1, h1); ii++) {
			pbuf[ii] = int( ii * f);
		}

		if (pbuf[h1 - 1] >= h) {
			pbuf[h1 - 1] = h - 1;
		}

		if (pbuf[w1 - 1] >= w) {
			pbuf[w1 - 1] = w - 1;
		}

		uchar *pimg = image.GetData();
		uchar *pdsimg = dsimage.GetData();
		for (ii = ind = 0; ii < h1; ii++) {
			i = pbuf[ii];

			for (jj = 0; jj < w1; jj++) {
				ind2 = (i * w + pbuf[jj]) * 3;
				pdsimg[ind++] = pimg[ind2++];
				pdsimg[ind++] = pimg[ind2++];
				pdsimg[ind++] = pimg[ind2++];
			}
		}

		delete []pbuf; pbuf = NULL;
	}
	else{
		resizeImageBilinear(image.GetData(), image.GetWidth(), image.GetHeight(), image.GetWidth()*3,
			dsimage.GetData(), dsimage.GetWidth(), dsimage.GetHeight(), dsimage.GetWidth()*3, 3);
	}
}

///////////////////////// image resize ///////////////////////////////
void resizeImageBilinear(const uchar* src, const int src_width, 
	const int src_height, const int src_step, uchar* dst, const int dst_width, 
	const int dst_height, const int dst_step, const int channels)
{
	// OpenCV stores the coefficients in short-type variables. 
	// each value is multiplied by INTER_RESIZE_COEF_SCALE (2048).
	const int INTER_RESIZE_COEF_BITS=11;
	const int INTER_RESIZE_COEF_SCALE=1 << INTER_RESIZE_COEF_BITS;

	const double scale_w = (double) src_width / (double) dst_width;
	const double scale_h = (double) src_height / (double) dst_height;

	int* src_w = new int[2 * dst_width];
	int* src_h = new int[2 * dst_height];
	short* coef_w = new short[2 * dst_width];
	short* coef_h = new short[2 * dst_height];

	for (int w = 0; w < dst_width; ++w)
	{
		float frac_w = (float)((w+0.5)*scale_w - 0.5);
		int integ_w = (int) std::floor(frac_w);
		frac_w -= integ_w;

		src_w[2*w] = clip(integ_w, 0, src_width);
		src_w[2*w+1] = clip(integ_w + 1, 0, src_width);
		coef_w[2*w] = saturate_cast_short((1.0f-frac_w)*INTER_RESIZE_COEF_SCALE);
		coef_w[2*w+1] = saturate_cast_short(frac_w*INTER_RESIZE_COEF_SCALE);
	}

	for (int h = 0; h < dst_height; ++h)
	{
		float frac_h = (float)((h+0.5)*scale_h - 0.5);
		int integ_h = (int) std::floor(frac_h);
		frac_h -= integ_h;

		src_h[2*h] = clip(integ_h, 0, src_height);
		src_h[2*h+1] = clip(integ_h + 1, 0, src_height);
		coef_h[2*h] = saturate_cast_short((1.0f-frac_h)*INTER_RESIZE_COEF_SCALE);
		coef_h[2*h+1] = saturate_cast_short(frac_h*INTER_RESIZE_COEF_SCALE);
	}

	for (int h = 0; h < dst_height; ++h)
	{
		for (int w = 0; w < dst_width; ++w)
		{
			for (int c = 0; c < channels; ++c)
			{
				int temp0 = coef_w[2*w] * src[src_h[2*h]*src_step + src_w[2*w]*channels + c] + 
					coef_w[2*w+1] * src[src_h[2*h]*src_step + src_w[2*w+1]*channels + c];
				int temp1 = coef_w[2*w] * src[src_h[2*h+1]*src_step + src_w[2*w]*channels + c] + 
					coef_w[2*w+1] * src[src_h[2*h+1]*src_step + src_w[2*w+1]*channels + c];
				dst[h*dst_step + w*channels + c] = 
					uchar(( ((coef_h[2*h] * ( temp0 >> 4)) >> 16) + ((coef_h[2*h+1] * (temp1  >> 4)) >> 16) + 2)>>2);
			}
		}
	}
	delete src_w;
	delete src_h;
	delete coef_w;
	delete coef_h;
}

// only support convention of .png
int convert(CImage &cImg, char* inname)
{
	static png_FILE_p fpin;
	static png_FILE_p fpout;

	if (strcmp(strrchr(inname,'.')+1, "png"))
		return -1;

	if(1){
		FILE* file = fopen(inname, "rb");

		png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);

		png_infop info_ptr = png_create_info_struct(png_ptr);

		setjmp(png_jmpbuf(png_ptr));

		png_init_io(png_ptr, file);

		png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_EXPAND, 0);

		int m_width = png_get_image_width(png_ptr, info_ptr);
		int m_height = png_get_image_height(png_ptr, info_ptr);

		int color_type = png_get_color_type(png_ptr, info_ptr);
		int pos = 0;
		png_bytep* row_pointers = png_get_rows(png_ptr, info_ptr);

		cImg.CreateImage(m_width, m_height, true);

		int ulChannels = png_get_channels(png_ptr, info_ptr);
		if(ulChannels==4 ){
			for(int i = 0; i < m_height; i++)
			{
				png_byte r, g, b, a;
				png_byte *dst = cImg.GetData() +  i * 3 * m_width;
				for(int j = 0; j < (ulChannels * m_width); j += ulChannels)
				{
					r = row_pointers[i][j + 2]; // blue
					g = row_pointers[i][j + 1]; // green
					b = row_pointers[i][j];   // red
					a = row_pointers[i][j + 3]; // alpha

					if(a==0){
						*dst++ = 255;
						*dst++ = 255;
						*dst++ = 255;
					}
					else{
						*dst++ = 255 - (255 - r) * a / 255; /* note the reverse order  */
						*dst++ = 255 - (255 - g) * a / 255;
						*dst++ = 255 - (255 - b) * a / 255;
					}
				}
			}
		}
		else if(ulChannels==3){
			for(int i = 0; i < m_height; i++)
			{
				png_byte r, g, b;
				png_byte *dst = cImg.GetData() +  i * 3 * m_width;
				for(int j = 0; j < (ulChannels * m_width); j += ulChannels)
				{
					r = row_pointers[i][j + 2]; // blue
					g = row_pointers[i][j + 1]; // green
					b = row_pointers[i][j];   // red

					*dst++ = r; /* note the reverse order  */
					*dst++ = g;
					*dst++ = b;
				}
			}
		}
		else if(ulChannels==2){
			for(int i = 0; i < m_height; i++)
			{
				png_byte r, a;
				png_byte *dst = cImg.GetData() +  i * 3 * m_width;
				for(int j = 0; j < (ulChannels * m_width); j += ulChannels)
				{
					a = row_pointers[i][j + 1]; // alpha
					r = row_pointers[i][j];   // red

					png_byte temp = 255 - (255 - r) * a / 255;
					*dst++ = temp;
					*dst++ = temp;
					*dst++ = temp;
				}
			}
		}
		else if(ulChannels==1){
			for(int i = 0; i < m_height; i++)
			{
				png_byte r, a;
				png_byte *dst = cImg.GetData() +  i * 3 * m_width;
				for(int j = 0; j < (ulChannels * m_width); j += ulChannels)
				{
					r = row_pointers[i][j];   // red

					png_byte temp = 255 - (255 - r) * a / 255;
					*dst++ = temp;
					*dst++ = temp;
					*dst++ = temp;
				}
			}
		}
		png_destroy_read_struct(&png_ptr, &info_ptr, 0);
		fclose(file);

	}
	else{

		png_structp read_ptr;
		png_infop read_info_ptr, end_info_ptr;

		png_bytep row_buf;
		png_uint_32 y;
		int num_pass, pass;
		png_uint_32 width, height;//宽度，高度
		int bit_depth, color_type;//位深，颜色类型
		int interlace_type, compression_type, filter_type;//扫描方式，压缩方式，滤波方式
		//读
		row_buf = NULL;
		//打开读文件
		if ((fpin = fopen(inname, "rb")) == NULL)
		{
			fprintf(stderr,"Could not find input file %s\n", inname);
			return (1);
		}

		//初始化1
		read_ptr =
			png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		read_info_ptr = png_create_info_struct(read_ptr);
		end_info_ptr = png_create_info_struct(read_ptr);
		//初始化2
		png_init_io(read_ptr, fpin);

		png_read_info(read_ptr, read_info_ptr);    
		//（1）IHDR
		png_get_IHDR(read_ptr, read_info_ptr, &width, &height, &bit_depth,&color_type, &interlace_type, &compression_type, &filter_type);

		//（2）cHRM
		//读取白色度信息  白/红/绿/蓝 点的x,y坐标，这里采用整形，不采用浮点数
		png_fixed_point white_x, white_y, red_x, red_y, green_x, green_y, blue_x,blue_y;

		png_get_cHRM_fixed(read_ptr, read_info_ptr, &white_x, &white_y,	&red_x, &red_y, &green_x, &green_y, &blue_x, &blue_y);
		//（3）gAMA
		png_fixed_point gamma;

		png_get_gAMA_fixed(read_ptr, read_info_ptr, &gamma);
		//（4）iCCP
		png_charp name;
		png_bytep profile;
		png_uint_32 proflen;

		png_get_iCCP(read_ptr, read_info_ptr, &name, &compression_type,	&profile, &proflen);

		//(5)sRGB
		int intent;
		png_get_sRGB(read_ptr, read_info_ptr, &intent);

		//(7)PLTE
		png_colorp palette;
		int num_palette;

		png_get_PLTE(read_ptr, read_info_ptr, &palette, &num_palette);

		//(8)bKGD
		png_color_16p background;

		png_get_bKGD(read_ptr, read_info_ptr, &background);

		//(9)hist

		png_uint_16p hist;

		png_get_hIST(read_ptr, read_info_ptr, &hist);

		//(10)oFFs
		png_int_32 offset_x, offset_y;
		int unit_type;

		png_get_oFFs(read_ptr, read_info_ptr, &offset_x, &offset_y,	&unit_type);

		//(11)pCAL
		png_charp purpose, units;
		png_charpp params;
		png_int_32 X0, X1;
		int type, nparams;

		png_get_pCAL(read_ptr, read_info_ptr, &purpose, &X0, &X1, &type, &nparams, &units, &params);
		//(12)pHYs

		png_uint_32 res_x, res_y;

		png_get_pHYs(read_ptr, read_info_ptr, &res_x, &res_y, &unit_type);
		//(13)sBIT
		png_color_8p sig_bit;

		png_get_sBIT(read_ptr, read_info_ptr, &sig_bit);

		//（14）sCAL
		int unit;
		png_charp scal_width, scal_height;

		png_get_sCAL_s(read_ptr, read_info_ptr, &unit, &scal_width,	&scal_height);

		//(15)iTXt
		png_textp text_ptr;
		int num_text;

		png_get_text(read_ptr, read_info_ptr, &text_ptr, &num_text);

		//(16)tIME,这里我们不支持RFC1123
		png_timep mod_time;

		png_get_tIME(read_ptr, read_info_ptr, &mod_time);

		//(17)tRNS
		png_bytep trans_alpha;
		int num_trans;
		png_color_16p trans_color;

		png_get_tRNS(read_ptr, read_info_ptr, &trans_alpha, &num_trans,	&trans_color);

		num_pass = 1;
		cImg.CreateImage(width, height, true);
		for (pass = 0; pass < num_pass; pass++)
		{
			for (y = 0; y < height; y++)
			{
				//分配内存
				row_buf = (png_bytep)png_malloc(read_ptr,
					png_get_rowbytes(read_ptr, read_info_ptr));
				png_read_rows(read_ptr, (png_bytepp)&row_buf, NULL, 1);

				//liu add
				png_byte *src, *dst;
				png_byte r, g, b;

				src = row_buf;
				dst = cImg.GetData() +  y * 3 * width;
				for (int xImg = 0; xImg < (int)width; xImg++)
				{
					r = *src++;
					g = *src++;
					b = *src++;		
					//a = *src++;

// 					if(a==0){
// 						*dst++ = 255;
// 						*dst++ = 255;
// 						*dst++ = 255;
// 					}
// 					else{
						*dst++ = b; //255 - (255 - b) * a / 255; /* note the reverse order  */
						*dst++ = g; //255 - (255 - g) * a / 255;
						*dst++ = r; //255 - (255 - r) * a / 255;
					//}
				}

				png_free(read_ptr, row_buf);
				row_buf = NULL;
			}
		}

		//结束
		//png_read_end(read_ptr, end_info_ptr);
		////
		////tTXt
		//png_get_text(read_ptr, end_info_ptr, &text_ptr, &num_text);

		////tIME
		//png_get_tIME(read_ptr, end_info_ptr, &mod_time);

		png_free(read_ptr, row_buf);
		row_buf = NULL;

		png_destroy_read_struct(&read_ptr, &read_info_ptr, &end_info_ptr);

		fclose(fpin);
	}

	return 1;
}