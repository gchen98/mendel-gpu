#include"../cl_constants.h"
#include"../io_manager.hpp"
#include"../mendel_gpu.hpp"
#include"../read_utilities.hpp"
#ifdef USE_GPU
#include"../cl_templates.hpp"
#endif

void DenovoReadsMendelGPU::init_opencl(){
  cerr<<"Initializing OpenCL for DenovoReadsMendelGPU\n";
  if(run_gpu){
#ifdef USE_GPU
    DenovoMendelGPU::init_opencl();
    read_penetrance->init_opencl();
  #endif
  }
}

void DenovoReadsMendelGPU::free_opencl(){
  if(run_gpu){
#ifdef USE_GPU
    read_penetrance->free_opencl();
    DenovoMendelGPU::free_opencl();
#endif
  }
}
