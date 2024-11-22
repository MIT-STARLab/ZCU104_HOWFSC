#include "appaccel.h"
#include <climits>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

std::vector<unsigned char> readFile(const char* fileName) 
{
    std::vector<unsigned char> buf;
    struct stat sb;
    if (stat(fileName, &sb) == -1) return buf; // file not found, return empty buffer on failure
    int fd = open(fileName, O_RDONLY);
	if (fd == -1) return buf; 
	
	printf("Loading %s\n",fileName);
	buf.resize(sb.st_size);
	read(fd, buf.data(),sb.st_size);
	close(fd);
    return buf;
}

bool listPlatforms(std::vector<cl::Platform> platforms)
{
	cl_int err;
    err = cl::Platform::get(&platforms);
	if (err != CL_SUCCESS) {printf("No platforms\n"); return false;}
    for (size_t j = 0; j < platforms.size(); j++) 
    {
        cl::Platform platform = platforms[j];
        std::string platformName = platform.getInfo<CL_PLATFORM_NAME>(&err);
        if (err == CL_SUCCESS) puts(platformName.c_str());
        puts("\n");
    }
    return true;
}

bool programXCL(const char* binaryFile, cl::Program & program,cl::Context & context,cl::CommandQueue & queue)
{
	cl_int err;
    std::vector<cl::Device> devices;
	if (!getDevices(devices, err)) return false;
	auto fileBuf = readFile(binaryFile);
	if (fileBuf.size())
	{
    	cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};

        for (unsigned int i = 0; i < devices.size(); i++) 
        {
            auto device = devices[i];
            // Creating Context and Command Queue for selected Device
            context = cl::Context(device, nullptr, nullptr, nullptr, &err);
            if (err != CL_SUCCESS) {printf("Error %d in cl::Context, device=%d\n",err,i); return false;}

            queue = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
            if (err != CL_SUCCESS) {printf("Error %d in cl::CommandQueue, device=%d\n",err,i); return false;}

            program = cl::Program(context, {device}, bins, nullptr, &err);
            if (err != CL_SUCCESS) 
            {
                printf("Error %d in cl::Program, device=%d\n",err,i);
            } 
            else 
            {
                printf("Device %d programed\n",i);;
                return true;
            }
        }
	}
    return false;
}

bool is_xpr_device(const char* device_name) 
{
    return strstr(device_name, "xpr") != nullptr;
}

bool is_emulation() 
{
    return (getenv("XCL_EMULATION_MODE") != nullptr);
}

bool is_hw_emulation() 
{
    char* xcl_mode = getenv("XCL_EMULATION_MODE");
    if (xcl_mode == nullptr) return false;
    if (!strcmp(xcl_mode, "hw_emu")) return true;
    return false;
}

bool getDevices(std::vector<cl::Device> &devices, cl_int &err, const char* vendorName) //defaults to Xilinx
{
    size_t i;
    std::vector<cl::Platform> platforms;
    err = cl::Platform::get(&platforms);
    if (err != CL_SUCCESS) {printf("No platforms\n"); return false;}
    cl::Platform platform;
    for (i = 0; i < platforms.size(); i++) 
    {
        platform = platforms[i];
        std::string platformName = platform.getInfo<CL_PLATFORM_NAME>(&err);
        if (err != CL_SUCCESS) {printf("Couldn't query platform name\n"); return false;}
        if (!(strcmp(platformName.c_str(),vendorName))) 
        {
            err = platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);
    		if(err == CL_SUCCESS){return true;};
        	printf("Error reading devices from %s, err %d\n",platformName.c_str(),err);
    		return false;
        }
    }
    return true;
}
