from personal_dna_analyzer.analyzer import main, get_arguments

if __name__ == '__main__':
    args = get_arguments()
    main(args.input, args.output, args.force_reload, args.data_dir)
