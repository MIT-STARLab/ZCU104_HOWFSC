
# FFT IP Core

## General Static Configuration Parameters Notes
The FFT IP Core supports both floating-point `float` and fixed-point `ap_fixed` data types. For optimal resource usage, it's recommended to use `ap_fixed` if you understand your input data requirements, as it consumes fewer resources.

### Input Width

- **Floating Point (`float`)**: Always 32-bit. Any other specified size is considered invalid.
- **Fixed Point (`ap_fixed`)**: Can be between 8 and 34 bits.

### Output Width
- **Floating Point (`float`)**: Always 32-bit. Any other specified size is considered invalid.
- **Fixed Point (`ap_fixed`)**: Depends on the scaling option. If scaled, then the same as the `input_width`. If not scaled, the output width is `input_width + max_nfft + 1`.



#### Notes
- Even though you can have different `port_width` values, if your kernel interacts with software define the variables with multiples of 8 bytes due to byte addressing constraints.
- For `ap_fixed`, scaled configurations must be normalized to `fix_x_[x-1]`, while unscaled configurations can be arbitrary.

#### Example

If `FFT_NFFT_MAX = log2(data_length) = 11` and `input_width = 32`, then:
- `output_width = input_width + FFT_NFFT_MAX + 1 = 43`

Define the data types as follows:
- `typedef ap_fixed<FFT_INPUT_WIDTH, 1> in_data_t = ap_fixed<32, 1>`
- `typedef ap_fixed<48, 17> out_data_t = ap_fixed<32 + (11 + 1) + 4 = 48, 1 + (11 + 1) + 4>`

Here, 4 extra bits are added to the width of the output data type to ensure it is a multiple of 8-bits.




### Max FFT Size (`max_nfft`)
The FFT data set size is specified as `1 << max_nfft`. `max_nfft` = `log2(Data length)`.



### Runtime Configurability (`has_nfft`)
The size can be configured at runtime, but must be a power of 2 and smaller than the maximum size. This flexibility requires additional resources.

### Architecture Optimization (`arch_opt`)
This parameter configures the implementation architecture, including pipeline and burst options. Higher I/O options increase throughput but trade-off with resource usage.
For more info check the AMD FFT IP Core document.

### Phase Factor Width (`phase_factor_width`)
Refers to the bit width used for representing twiddle factors. Configure this for internal phase factor precision.

### Output Ordering (`ordering_opt`)
Default is reversed (bit-reversed). Reversing bits affects burst time and pipeline RAM usage.

### Scaling Options (`scaling_opt`)

- **[Scaling Schedule Supporting File](https://support.xilinx.com/s/article/1160838?language=en_US)**
  
  FFT at each stage can increase numbers by a couple of bits, so the output might need extra bits to avoid overflow:
  
  - **Option 1**: Set the output data type width to `input_width + max_nfft + 1` to avoid scaling issues.
  - **Option 2**: Provide a scaling schedule to scale data at each FFT stage. Note:
    - The scaling schedule is passed at runtime.
    - Too much early scaling can lose precision; too little can cause overflow.
    - Conservative scaling schedules are provided in [PG109].
      - For N=1024
        - Radix-4 Burst I/O  :  [10 10 10 10 11]
        - Pipelined Stream I/O : [10 10 10 10 11]  For N=2048 : [01 10 10 10 10 11]
        - Radix-2 Burst I/O or Lie : [01 01 01 01 01 01 01 01 01 10]
  - **Option 3**: Let the tool automatically determine scaling, though this may use more resources.

  For floating-point data, scaling is handled internally and a scaling schedule is not required (SCALE_SCH is ignored).

  

### More Notes
- Floating-Point format is not available in multichannel configurations.


