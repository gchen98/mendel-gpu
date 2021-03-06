#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void DenovoMendelGPU::impute_genotypes_opencl(){
  int c_snp = g_center_snp_start-g_left_marker;
  int current_snp = g_center_snp_start;
  cerr<<"begin diploid impute_genotypes at SNP "<<current_snp<<", center: "<<c_snp<<endl;
  //return ;
  bool debug_dosage = false;
  bool debug_geno = false;
  bool debug_posterior = g_left_marker<-1;
  bool debug_pen = false;
  int penmat_person = 0;
  int max_geno=g_genotype_imputation?3:4;
  cerr<<"Max geno "<<max_geno<<endl;
  for(int i=0;i<g_max_haplotypes;++i){
    center_dosage[i] = g_haplotype[i*g_max_window+c_snp];
  }
  float subject_dosages[g_people];
  if (run_gpu){
    float subject_posterior[4*g_people];
    int subject_genotypes[g_people];
    double start = clock();
    #ifdef USE_GPU
    writeToBuffer(buffer_center_dosage, g_max_haplotypes, center_dosage, "buffer_center_dosage");
    writeToBuffer(buffer_hap_perm, g_max_haplotypes, g_hap_perm, "buffer_hap_perm");
    runKernel("kernel_impute_genotype_denovo",kernel_impute_genotype_denovo,g_people*BLOCK_WIDTH,1,1,BLOCK_WIDTH,1,1);
    cerr<<"Launched impute_genotype_denovo\n";
    readFromBuffer(buffer_subject_posterior_prob, g_people*4,subject_posterior,"buffer_subject_posterior_prob");
    float posterior[g_people*max_geno];
    for(int i=0;i<g_people;++i){
      for(int j=0;j<max_geno;++j){
        posterior[i*max_geno+j] = subject_posterior[i*4+j];
      }
    }
    if (debug_pen){
      float * debug_pen = new float[g_people*penetrance_matrix_size];
      readFromBuffer(buffer_penetrance_cache,g_people*penetrance_matrix_size,debug_pen,"buffer_penetrance_cache");
      for(int i=0;i<g_max_haplotypes;++i){
        for(int j=0;j<g_max_haplotypes;++j){
if (g_active_haplotype[i] && g_active_haplotype[j]){
          cerr<<i<<","<<j<<":"<<debug_pen[penmat_person*penetrance_matrix_size+i*g_max_haplotypes+j]<<endl;
}
        }
      }
      delete[]debug_pen;
    }
    if (debug_posterior){
      for(int i=0;i<g_people;++i){
        float denom = 0;
        for(int j=0;j<max_geno;++j){
          denom+=subject_posterior[i*4+j];
        }
        cerr<<"GPU_POSTERIOR "<<i<<","<<current_snp;
        for(int j=0;j<max_geno;++j){
          cerr<<" "<< subject_posterior[i*4+j];
        }
        cerr<<endl;
      }
    }
    readFromBuffer(buffer_subject_genotype, g_people,subject_genotypes,"buffer_subject_genotype");
    if (debug_geno){
      for(int i=0;i<g_people;++i){
        cerr<<"GPU GENO:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_genotypes[i]<<endl;
      }
    }
    if (g_genotype_imputation){
      readFromBuffer(buffer_subject_dosage, g_people,subject_dosages,"buffer_subject_dosage");
    }
    if (debug_dosage){
      for(int i=0;i<g_people;++i){
        cerr<<"GPU_DOSAGE:\t"<<i<<"\t"<<current_snp<<"\t"<<subject_dosages[i]<<endl;
      }
    }
    if (max_geno==3){
      io_manager->writeDosage(current_snp,subject_dosages,g_people);
      float rsq = compute_rsq(subject_dosages,1,0);
      io_manager->writeQuality(current_snp,rsq);
    }
    for(int i=0;i<g_people;++i){
      int hap1_allele = true_haps[(2*i)*g_snps+current_snp];
      int hap2_allele = true_haps[(2*i+1)*g_snps+current_snp];
      subject_genotypes[i] = hap1_allele*2+hap2_allele;
    }
    io_manager->writeGenotype(current_snp,subject_genotypes,g_people);
    io_manager->writePosterior(max_geno,current_snp,posterior,g_people);
    cerr<<"Elapsed time: "<<(clock()-start)/CLOCKS_PER_SEC<<endl;
    #endif
  }
}
