#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <cstring>

#include "helper.h"
#include "bgen_parser.h"

using std::string;
using namespace std;

int main(int argc, char* argv[])
{
	bool VERBOSE=0;
	ifstream bgenfile("/lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/001/Stripe/ukb_imp_chr22_v2.bgen", ios::in | ios::binary);
	file_info bgen_finfo;
	bgen_finfo=parse_header(bgenfile, 0);
	ofstream out("/nfs/users/nfs_a/ag15/lol", ios :: out | ios :: binary);
	write_header(out, bgen_finfo);
	bgenfile.close();
	exit(0);
	unsigned short id_len;
	string id;
	string rsid;
	string chr;
	uint32_t pos;

	while(!bgenfile.eof()){
		bgenfile.read(reinterpret_cast<char *>(&id_len), sizeof(id_len));
		if(VERBOSE){info("First variant ID length: ", id_len);}
		char varid[id_len];
		bgenfile.read(varid, id_len);
		id=std::string(varid);
		if(VERBOSE){info("Variant ID: ", id);}

		bgenfile.read(reinterpret_cast<char *>(&id_len), sizeof(id_len));
		if(VERBOSE){info("First variant rsID length: ", id_len);}
		char varrs[id_len];
		bgenfile.read(varrs, id_len);
		rsid=std::string(varrs);
		if(VERBOSE){info("Variant rsID: ", rsid);}

		bgenfile.read(reinterpret_cast<char *>(&id_len), sizeof(id_len));
		if(VERBOSE){info("First variant chr length: ", id_len);}
		char varchr[id_len];
		bgenfile.read(varchr, id_len);
		chr=std::string(varchr);
		if(VERBOSE){info("Variant chr: ", chr);}

		bgenfile.read(reinterpret_cast<char *>(&pos), sizeof(pos));
		if(VERBOSE){info("Variant position: ", pos);}

	//alleles
		unsigned short num_alleles;
		bgenfile.read(reinterpret_cast<char *>(&num_alleles), sizeof(num_alleles));
		if(VERBOSE){info("Variant has ", num_alleles, "alleles");}
		string alleles[num_alleles];
		for(int i=0;i<num_alleles;i++){
			uint32_t len_allele;
			bgenfile.read(reinterpret_cast<char *>(&len_allele), sizeof(len_allele));
			char allele_c[len_allele];
			bgenfile.read(allele_c, len_allele);
			string allele;
			alleles[i] = std::string(allele_c);
			if(VERBOSE){info("\t allele ", allele, "(size ", len_allele, ")");}
		}

	//genotypes
	uint32_t len_gblock; // compressed size
	uint32_t len_ugblock; // uncompressed size
	bgenfile.read(reinterpret_cast<char *>(&len_gblock), sizeof(len_gblock));
	bgenfile.read(reinterpret_cast<char *>(&len_ugblock), sizeof(len_ugblock));
	long unsigned int len_ugblock_l=(long unsigned int)len_ugblock;
	long unsigned int len_gblock_l=(long unsigned int)len_gblock;

	char gblock[len_gblock-4]; // read compressed block
	bgenfile.read(gblock, len_gblock-4);
	
	uint8_t dest[len_ugblock];
	int unzip_code;

	unsigned char* gblock_u=reinterpret_cast<unsigned char*>(gblock);
	if(VERBOSE){info("Compressed genotype block of size: ", len_gblock_l);}
	if(VERBOSE){info("Uncompressed genotype block  size: ", len_ugblock_l);}
	unzip_code=uncompress(dest, &len_ugblock_l, gblock_u, len_gblock_l);
	if(VERBOSE){info("Decompress returned ",unzip_code);}

	uint8_t* temp=&dest[0];
	uint8_t snum[4];
	memcpy(&snum, temp, 4);
	uint32_t ns=*reinterpret_cast<uint32_t*>(snum);
	if(VERBOSE){info("Number of samples in this variant: ",ns);}

	temp+=4;
	uint8_t numall_p[2];
	memcpy(numall_p, temp, 2);
	int numall=*reinterpret_cast<uint16_t*>(numall_p);
	if(VERBOSE){info("Number of alleles for this variant: ",numall);}

	temp+=2;
	uint8_t maxpl[1];
	memcpy(maxpl, temp, 1);
	if(VERBOSE){info("Min ploidy for this variant: ",(int)*maxpl);}

	temp+=1;
	uint8_t minpl[1];
	memcpy(minpl, temp, 1);
	if(VERBOSE){info("Max ploidy for this variant: ",(int)*minpl);}

	// reading ploidy and missingness information
	temp+=1;
	uint8_t ploidy[ns];
	memcpy(ploidy, temp, ns);

	// reading phasedness
	temp+=ns;
	uint8_t phased[1];
	memcpy(phased, temp, 1);
	if(VERBOSE){info("Data is phased: ",(int)*phased);}

	// reading number of bits per genotype
	temp+=1;
	uint8_t nbits[1];
	memcpy(nbits, temp, 1);
	if(VERBOSE){info("Data is stored in  ",(int)*nbits, "bits.");}

	float flnb=(float) (2<< nbits[0] -1);

	temp+=1;

	// probabilities are stored depending on ploidy and allele number per sample.
	// they cannot be read in one go.
	
	// genotype assignment loop
	bool missingness[ns];
	unsigned long int totmiss=0;
	for(int i = 0; i < ns; ++i){
		uint8_t prob;
		missingness[i]=ploidy[i] & 128;
		totmiss+=missingness[i];
		ploidy[i] <<= 1;
		ploidy[i] >>= 1;
		int numprob=(ploidy[i]+numall-1)*choose(numall-1, numall-1)-1;

		// read genotype
		// the following assumes that probabilities are always stored in 8-bit rep, i.e. nbits==8 always
		uint8_t probs[numprob];
		memcpy(probs, temp, numprob);
		temp+=numprob;
		float sump=1;
		float dosage=0;
		for (unsigned int j=0; j<numprob; j++){
			float flpr=(float) probs[j];
			flpr=flpr/flnb;
			sump-=flpr;
			dosage+=(float)j*flpr;
			//cout << flpr<< ":";
		}
		dosage+=2*sump;
		//cout << dosage <<" ";

	}
	//cout << "\n";
	if(VERBOSE){info("Total ploidy is ", ploidy);}


	// float fns=(float)ns*2;
	// float fpl=(float)ploidy;
	// float fms=(float)totmiss;
	// std::cout << id<<"\t"<<rsid<<"\t" << (int)*minpl <<"\t" <<fns<<"\t"<< fpl<<"\t"<< fms <<"\t"<<fpl/(fns-totmiss) <<"\n";
}
}

	// char zipped_block_c[len_gblock];
	// bgenfile.read(zipped_block_c, len_gblock);
	// string zipped_block;
	// zipped_block=std::string(zipped_block_c);
	// info(uncompress(zipped_block_c, 0));
