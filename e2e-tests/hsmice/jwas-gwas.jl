using CSV, DataFrames, JWAS
using Pipe: @pipe

if length(ARGS) != 3
    print(stderr, "Usage: julia jwas-gwas.jl GENOTYPE-FILE PHENOTYPE-FILE OUTPUT-GWAS-FILE\n")
    exit(1)
end
genotype_file, phenotype_file, output_gwas_file = ARGS

genotype = CSV.read(genotype_file, DataFrame; drop=[:reference])
phenotype = @pipe(CSV.read(phenotype_file, DataFrame)
                  |> select(_,
                            "sample-id" => :sample_id,
                            "Biochem.ALP.resid" => :alp))

genotypes = get_genotypes(@pipe(
    hcat(DataFrame(snp_id=string.(genotype.chromosome, ":", genotype.position)),
         @pipe(genotype
               |> select(_, Not([:chromosome, :position]))))
    |> permutedims(_, :snp_id)
    |> rename(_, :snp_id => :sample_id)),
                          1.8e-5,
                          G_is_marker_variance=true,
                          quality_control=false,
                          center=false,
                          Pi=0.999,
                          estimateScale=true)

model = build_model("alp = intercept + genotypes")
out = runMCMC(model,
              phenotype,
              Pi=0.99,
              estimateScale=true,
              chain_length=50000,
              burnin=5000,
              output_heritability=true,
              output_samples_frequency=100,
              # fix seed for reproducibility
              seed=1)

gwas = GWAS("results/MCMC_samples_marker_effects_genotypes_alp.txt")
CSV.write(output_gwas_file,
          transform(gwas,
                    :marker_ID => ByRow(x -> split(x, ":")[1]) => :chromosome,
                    :marker_ID => ByRow(x -> parse(Int, split(x, ":")[2])) => :position);
          delim="\t")
