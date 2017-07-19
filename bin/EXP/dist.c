#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "errno.h"
#include "pthread.h"

//XXX 计数节点，该节点主要用于feature的计数，每个feature只有一个count_node节点 XXX//
typedef struct _count_node_
{
	char mutex;		//多线程互斥锁，暂时没有考虑原子操作
	char strand;	//feature方向
	int key;		//保存最近统计的read的起始结束位置之和
	int feature;	//feature类型
	int left;		//feature左坐标
	int right;		//feature右坐标
	int reads_num;	//落在该feature的reads数量的2倍
	int cov;		//在该feature的所有reads的碱基数
	struct _count_node_ * next;	//链表的下一个节点
}count_node;

//XXX 索引节点，该节点主要用于快速检索一个索引里的所有feature的位置信息 XXX//
typedef struct _node_
{
	char strand;
	int feature;
	int left;
	int right;
	count_node * count;		//该index中的feature指向count_node计数节点的指针
	struct _node_ * next;
}node;

typedef struct _arg_
{
	int step;
	int map_mode;
	long start;
	long end;
	node *** search_list;
	char * sam_file;
}arg;

void usage(char * program_name)
{
	printf("usage:%s sam_file chr_file anno_file out_file step_length [thread_num]\n",program_name);
	exit(0);
	return;
}

int cal_chrom(int * chr_len,char * chr_file)
{
	int chr_num = 0;
	FILE * fp = fopen(chr_file,"r");
	if(fp == NULL)
	{
		printf("open %s error\n",chr_file);
		exit(1);
	}
	char * buf = (char *)malloc(1024 * sizeof(char));

	char * chr = (char *)malloc(8 * sizeof(char));
	int length = 0;
	while(fgets(buf,1024,fp))
	{
//Chr1	TAIR10	chromosome	1	30427671	.	.	.	ID=Chr1;Name=Chr1
		if( 2 < sscanf(buf,"%s\t%d",chr,&length))
		{
			continue;
		}
		*(chr_len + *((unsigned short *)(chr + 3))) = length;
	}
	free(buf);
	free(chr);
	fclose(fp);
	return chr_num;
}

///////////////////////////////////////////////////////////
//						建立查找表						 //
///////////////////////////////////////////////////////////
int build_search_list(char * anno_file,node *** search_list,int step,count_node ** count_list)
{
	//XXX 打开注释信息文件 XXX//
	FILE * fp = fopen(anno_file,"r");
	if(fp == NULL)
	{
		printf("open %s error\n",anno_file);
		exit(1);
	}
	char * buf = (char *)malloc(1024 * sizeof(char));		//行缓存

	char * chr = (char *)malloc(8 * sizeof(char));			//
	char * gene_id = (char *)malloc(128 * sizeof(char));	//
	char strand;					//方向 +/-
	long feature_count[64] = {0};	//用于统计每一个feature出现频率的表
	int feature = 0;				//注释文件feature的编码
	int start = 0;					//注释文件起始点
	int end = 0;					//注释文件结束点
	int end_point = 0;				//将feature添加进查找表的循环结束位点
	int i = 0;						//循环的index

	count_node * count_list_end[65536] = {0};			//每个feature的count节点的链表末尾
	bzero(count_list_end,65536 * sizeof(count_node *));
	//XXX 初始化feature的count节点链表 XXX//
	for(i = 0;i < 65536;i++)
	{
		*(count_list + i) = (count_node *)malloc(sizeof(count_node));
		(*(count_list + i))->next = NULL;
		count_list_end[i] = *(count_list + i);
	}

	node * new_node_ptr = NULL;		//新建node类型节点的指针变量
	count_node * count_ptr = NULL;	//新建count_node类型节点的指针变量
	while(fgets(buf,1024,fp))
	{
//Chr1	TAIR10	chromosome	1	30427671	.	.	.	ID=Chr1;Name=Chr1
		if( 6 < sscanf(buf,"%s\t%d\t%c\t%d\t%d\t%s",chr,&feature,&strand,&start,&end,gene_id))
		{
			printf("sscanf error\n");
			exit(1);
		}
		feature_count[feature]++;		//统计feature出现的种类
		//XXX 为该feature创建一个计数的节点并初始化 XXX//
		count_ptr = (count_node *)malloc(sizeof(count_node));
		count_ptr->left = start;
		count_ptr->right = end;
		count_ptr->strand = strand;
		count_ptr->feature = feature;
		count_ptr->reads_num = 0;
		count_ptr->cov = 0;
		count_ptr->next = NULL;
		count_list_end[*((unsigned short *)(chr + 3))]->next = count_ptr;
		count_list_end[*((unsigned short *)(chr + 3))] = count_ptr;
		end_point = (end / step) * step;
		for(i = start;i < end_point;i += step)
		{
			//XXX 该feature在索引中的检索节点 XXX//
			new_node_ptr = (node *)malloc(sizeof(node));
			new_node_ptr->left = start;
			new_node_ptr->strand = strand;
			new_node_ptr->feature = feature;
			new_node_ptr->right = end;
			new_node_ptr->count = count_ptr;
			new_node_ptr->next = *(*(search_list + *((unsigned short *)(chr + 3))) + (i / step));
			*(*(search_list + *((unsigned short *)(chr + 3))) + (i / step)) = new_node_ptr;
//			printf("%d\t%d\t%p\n",*((unsigned short *)(chr + 3)),(i / step),new_node_ptr);
		}
		new_node_ptr = (node *)malloc(sizeof(node));
		new_node_ptr->left = start;
		new_node_ptr->strand = strand;
		new_node_ptr->feature = feature;
		new_node_ptr->right = end;
		new_node_ptr->count = count_ptr;
		new_node_ptr->next = *(*(search_list + *((unsigned short *)(chr + 3))) + (end / step));
		*(*(search_list + *((unsigned short *)(chr + 3))) + (end / step)) = new_node_ptr;
//		printf("%d\t%d\t%p\n",*((unsigned short *)(chr + 3)),(i / step),new_node_ptr);
	}
	int max_feature_num = 0;
	for(i = 0;i < 64;i++)
	{
		if(feature_count[i] > 0)
		{
			max_feature_num = i;
		}
	}
	fclose(fp);
	free(buf);
	free(gene_id);
	free(chr);
	return max_feature_num;
}

///////////////////////////////////////////////////////////
//			将每个read/fragment添加到查找表里			 //
///////////////////////////////////////////////////////////
inline void push_search_list(unsigned short chr,int f,int r,node *** search_list,int map_mode,char strand,int step,int reads_num)
{	//chr是read所在染色体的短整型编码，f是read起始，r是结束，map_mode在这里是全部read或者fragment//
	node * pointer = NULL;
	count_node * count = NULL;
	if((r / step) > (f / step))		//这个read/fragment跨过了一个window
	{
		pointer = *(*(search_list + chr) + (f / step));
		while(pointer != NULL)
		{
			count = pointer->count;
			while(count->mutex) free(0);
			count->mutex = 1;
			//XXX read完全包含在该feature里 XXX//
			if(pointer->left <= f && pointer->right >= r && count->key != reads_num)
			{
				if(map_mode)
				{
					count->reads_num += 1;
				}
				else
				{
					count->reads_num += 2;
				}
				count->cov += r - f + 1;
			}
			count->key = reads_num;
			count->mutex = 0;
			pointer = pointer->next;
		}
	}
	pointer = *(*(search_list + chr) + (r / step));
	while(pointer != NULL)
	{
		count = pointer->count;
		while(count->mutex) free(0);
		count->mutex = 1;
		if(pointer->left <= f && pointer->right >= r && count->key != reads_num)
		{
			if(map_mode)
			{
				count->reads_num += 1;
			}
			else
			{
				count->reads_num += 2;
			}
			count->cov += r - f + 1;
		}
		count->key = reads_num;
		count->mutex = 0;
		pointer = pointer->next;
	}
	return;
}

///////////////////////////////////////////////////////////
//			线程函数统计每个feature的reads表达			 //
///////////////////////////////////////////////////////////
void * statistics_pthread(void * argument)
{
	//XXX 初始化参数 XXX//
	arg * thread_arg = (arg *)argument;
	long start = thread_arg->start;
	long end = thread_arg->end;
	int step = thread_arg->step;
	int map_mode = thread_arg->map_mode;
	char * sam_file = thread_arg->sam_file;
	node *** search_list = thread_arg->search_list;
	//printf("start:%ld\tend:%ld\t%p\n",start,end,search_list);
	
	//XXX 打开文件定位文件指针到该线程处理的位置 XXX//
	FILE * fp = fopen(sam_file,"r");
	if(fp == NULL)
	{
		printf("open %s error\n",sam_file);
		exit(1);
	}
	fseek(fp,start,SEEK_SET);

	//XXX 初始化循环相关变量 XXX//
	unsigned short chr = 0;			//染色体字符串最后两个字符强制转换成短整型的变量
	unsigned int flag = 0;			//sam文件标记位
	int reads_num = 0;				//reads序号
	int pos = 0;					//sam文件中read的染色体位置
	int f,r;						//添加到查找表里的read片段前后端点
	int other;						//sscanf中无用的变量
	int left,right,gap = 0;			//sam文件中fragment mapping的左右两条reads起始位点和之间的间隔
	char strand;					//方向
	long nega = 0;					//临时变量
	long posi = 0;					//临时变量
	char * ptr = NULL;				//字符串指针
	char * buf = (char *)malloc(1024 * sizeof(char));		//sam文件的每一行的缓存
	char * map_info = (char *)malloc(32 * sizeof(char));	//map信息的字符串信息
	char * Chr = (char *)malloc(16 * sizeof(char));			//染色体字符串信息
	while(fgets(buf,1024,fp))
	{
		reads_num++;
		ptr = buf;
		if(ftell(fp) > end)
		{
			break;
		}
		while(*ptr != '\t') ptr++;ptr++;	//跳过read名
// SBS_0033_FC70L9TAAXX:1:4:16782:18489#0[2307316] 0       Chr5    1839926 255     28M182N52M      *       0       0
		if(5 < sscanf(ptr,"%d\t%s\t%d\t%d\t%s\t",&flag,Chr,&pos,&other,map_info)) exit(1);
		if(flag & 16)		//判断方向
		{
			strand = '-';
			nega++;
		}
		else
		{
			strand = '+';
			posi++;
		}
		chr = *((unsigned short *)(Chr + 3));
		//XXX 如果是fragment mapped XXX//
		if(strlen(map_info) > 5)
		{
			sscanf(map_info,"%dM%dN%dM",&left,&gap,&right);
			f = pos;
			r = pos + left - 1;
			push_search_list(chr,f,r,search_list,1,strand,step,reads_num);	//将左边的fragment添加到表里
			f = r + gap + 1;
			r = f + right - 1;
			push_search_list(chr,f,r,search_list,1,strand,step,reads_num);	//将右边的fragment添加到表里
		}
		else	//非fragment mapped
		{
			sscanf(map_info,"%dM",&left);
			f = pos;
			r = pos + left - 1;
			push_search_list(chr,f,r,search_list,0,strand,step,reads_num);	//将整个read添加到表里
		}
	}
	fclose(fp);
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////
//                               计算每个feature的表达                        //
////////////////////////////////////////////////////////////////////////////////
int cal_expr(node *** search_list,char * sam_file,int step,int threads_num,int map_mode)
{
	FILE * fp = fopen(sam_file,"r");
	if(fp == NULL)
	{
		printf("open %s error\n",sam_file);
		exit(1);
	}
	char * buf = (char *)malloc(1024 * sizeof(char));					//行缓存
	//XXX 分割sam文件以便多线程处理 XXX//
	fseek(fp,0,SEEK_END);
	long sam_file_size = ftell(fp);										//sam文件的大小
	long part_size = sam_file_size / threads_num;						//分成多份之后每一份的大小
	long * part_pos = (long *)malloc((threads_num + 1) * sizeof(long));	//精确分割位置的数组指针
	bzero(part_pos,threads_num * sizeof(long));
	int i = 0;

	//XXX 去除sam文件开始部分以@起始的行 XXX//
	long pos = 0;
	fseek(fp,0,SEEK_SET);
	while(fgets(buf,1024,fp))
	{
		if(*buf != '@')
		{
			break;
		}
		pos = ftell(fp);
	}

	//XXX 确定文件分割的精确位置 XXX//
	for(i = 0;i < threads_num - 1;i++)
	{
		*(part_pos + i) = pos;
		pos += part_size;
		fseek(fp,pos,SEEK_SET);
		fgets(buf,1024,fp);		//读完该行剩下的字符，以便下行从行首开始
		pos = ftell(fp);
	}
	*(part_pos + i) = pos;
	*(part_pos + i + 1) = sam_file_size;
	//printf("pos:%ld\tsam_file_size:%ld\n",pos,sam_file_size);

	pthread_t * thread_handle = (pthread_t *)malloc(threads_num * sizeof(pthread_t));
	arg * thread_args = (arg *)malloc(threads_num * sizeof(arg));
	bzero(thread_handle,threads_num * sizeof(pthread_t));
	bzero(thread_args,threads_num * sizeof(arg));

	//XXX 启动多线程处理 XXX//
	for(i = 0;i < threads_num;i++)
	{
		//为线程函数的参数结构体赋值
		(thread_args + i)->start = *(part_pos + i);
		(thread_args + i)->end = *(part_pos + i + 1);
		(thread_args + i)->sam_file = sam_file;
		(thread_args + i)->step = step;
		(thread_args + i)->map_mode = map_mode;
		(thread_args + i)->search_list = search_list;
		//启动一个线程
		pthread_create(thread_handle + i,NULL,statistics_pthread,thread_args + i);
	}
	//XXX 回收线程 XXX//
	for(i = 0;i < threads_num;i++)
	{
		pthread_join(*(thread_handle + i),NULL);
	}

	free(thread_args);
	free(thread_handle);
	fclose(fp);
	return 0;
}

void output_result(node *** search_list,int max_feature_num,int * chr_len,int step,count_node ** count_list,char * out_root)
{
	//XXX 初始化变量 XXX//
	unsigned short i = 0;
	count_node * ptr = NULL;	//计数节点的指针
	char c1,c2 = '\0';			//将染色体编码短整型翻译成两个字符型所用的变量
	char chr[8] = {0};			//染色体字符数组
	char filename[128] = {0};	//文件名数组，用于保存每个feature

	//XXX 创建多个文件，每个文件为对一个feature的统计 XXX//
	FILE ** file_ptr = (FILE **)malloc(max_feature_num * sizeof(FILE *));
	bzero(file_ptr,max_feature_num * sizeof(FILE *));
	for(i = 0;i < max_feature_num + 1;i++)
	{
		sprintf(filename,"%s.%d",out_root,i);
		*(file_ptr + i) = fopen(filename,"w");
		if(*(file_ptr + i) == NULL)
		{
			printf("create file %s error\n",filename);
		}
	}

	//XXX 对每个染色体进行输出，染色体编码方式在本文件中有说明 XXX//
	for(i = 0;i < 65535;i++)
	{
		c1 = (char)(i >> 8);
		c2 = (char)(i & 0xff);
		sprintf(chr,"chr%c%c",c2,c1);
		ptr = (*(count_list + i))->next;
		if(*(chr_len + i) > 0)
		{
			while(ptr)
			{
				//printf("flag\n");
				fprintf(*(file_ptr + ptr->feature),"%s\t%d\t%c\t%d\t%d\t%d\t%d\n",chr,ptr->feature,ptr->strand,ptr->left,ptr->right,ptr->reads_num,ptr->cov);
				ptr = ptr->next;
			}
		}
	}
	for(i = 0;i < max_feature_num;i++)
	{
		fclose(*(file_ptr + i));
	}
	return;
}

int main(int argc,char ** argv)
{
	//XXX 初始化参数 XXX//
	if(argc < 5)
	{
		usage(argv[0]);
	}
	int step = atoi(argv[5]);	//index的步长
	int threads_num = 2;		//线程数量
	char * anno_file = (char *)malloc(256 * sizeof(char));		//注释信息文件
	char * sam_file = (char *)malloc(256 * sizeof(char));		//sam文件
	char * chr_file = (char *)malloc(256 * sizeof(char));		//染色体长度文件
	char * out_file = (char *)malloc(256 * sizeof(char));		//输出文件，暂未用到
	strcpy(sam_file,argv[1]);
	strcpy(chr_file,argv[2]);
	strcpy(anno_file,argv[3]);
	strcpy(out_file,argv[4]);
	if(argv[6] != NULL)
	{
		threads_num = atoi(argv[6]);
	}

	int map_mode = 0;		//比对模式，暂时不使用
	//XXX 65536指的是一个无符号短整型变量的取值范围，同时也是两个字符型变量取值范围 XXX//
	int * chr_len = (int *)malloc(65536 * sizeof(int));								//染色体长度数组
	node *** search_list = (node ***)malloc(65536 * sizeof(node **));				//二维指针变量查找表，线性加链接
	count_node ** count_list = (count_node **)malloc(65536 * sizeof(count_node *));	//二维count节点链表
	bzero(chr_len,65536 * sizeof(int));
	bzero(search_list,65536 * sizeof(node **));
	bzero(count_list,65536 * sizeof(count_node *));

	printf("start:");
	fflush(stdout);
	system("date");
	cal_chrom(chr_len,chr_file);	//获取染色体长度
	int i = 0;
	for(i = 0;i < 65536;i++)		//为每条染色体查找表分配空间
	{
		if(*(chr_len + i))
		{
			*(search_list + i) = (node **)malloc((*(chr_len + i) / step) * sizeof(node *));
			bzero(*(search_list + i),(*(chr_len + i) / step) * sizeof(node *));
			//printf("%d\n",*(chr_len + i));
		}
	}

	int max_feature_num = build_search_list(anno_file,search_list,step,count_list);	//建立查找表
	/*//////////////debug//////////////////
	int j = 0;
	node * ptr = NULL;
	for(i = 0;i < 65536;i++)
	{
		if(*(chr_len + i))
		{
			for(j = 0;j < *(chr_len + i);j += step)
			{
				ptr = *(*(search_list + i) + (j / step));
				if(ptr) printf("<%d:%d>",i,j / step);
				while(ptr)
				{
					printf("->%d:%d",ptr->left,ptr->right);
					ptr = ptr->next;
				}
				if(*(*(search_list + i) + (j / step))) printf("\n");
			}
		}
	}
	*////////////////debug/////////////////
	cal_expr(search_list,sam_file,step,threads_num,map_mode);	//计算每个feature的表达量

	printf("output result:");
	fflush(stdout);
	system("date");
	output_result(search_list,max_feature_num,chr_len,step,count_list,out_file);		//输出结果
	printf("end:");
	fflush(stdout);
	system("date");
	return 0;
}
