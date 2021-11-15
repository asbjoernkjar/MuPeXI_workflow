

#this function splits a vcg file --> does only produce 1 of the splits 
#if n= 2 and splits = 100 --> will cut out mut 101-200

def split_vcf(in_path,out_path, start, end):
    #start inclusive
    #end non inclusive

    big_vcf = open(in_path, "r")
    lines = big_vcf.readlines()
    big_vcf.close()
    out = open(out_path, "w")
    start, end = int(start), int(end)
    
    c = 1
    print("splitting to sample:", out_path)
    print("interval:" , start, end)
    

    for line in lines:
        if line[0] == "#":
            out.write(line)
            continue
        line_split = line.strip().split("\t")
        if line_split[6] == "PASS":
            if c > end:
                break

            elif c > start:
                out.write(line)
            
            c += 1
    print("ended at:", c)
    print()
    out.close()


if __name__ == "__main__":
    split_vcf(snakemake.input.inp, snakemake.output.out, snakemake.params.start, snakemake.params.end)


