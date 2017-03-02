# -*- coding: utf-8 -*-
require "ffc/version"
require "ffc/cli"

require 'mkmf'
require 'systemu'
require "open3"
require "config"

module Conv

  class Conv
    attr_accessor :conv_items, :package_manager

    def initialize
      path_name = Pathname.new( File.expand_path(File.join(File.dirname(__FILE__), 'db.bin')))
      if !path_name.exist?
        self.save!
      end
      path_name.open('rb') do |f|
        @conv_items = Marshal.load(f)
        @conv_items.default_proc = ->(h, k) { h[k] = [] }
      end
    end

    def self.save!
      @conv_items = Hash.new{ |h,k| h[k] = [] }
      #items.each{ |item| item.input.each { |t| @conv_items[t] << item }}
      self.rehash.each{ |item| @conv_items[item.input] << item }
      @conv_items.default_proc = nil
      File.open(File.expand_path(File.join(File.dirname(__FILE__), 'db.bin')), 'w') {|f| f.write(Marshal.dump(@conv_items)) }
    end

    def self.update!(new_item)
      #TODO(NOT IMPLEMENTED)
      #new_item.last_access = Time.now
      #p new_item
      #p @conv_items.select{ |item| item.command == new_item.command }#[0] = new_item
      #@conv_items.default_proc = nil
      #File.open(File.expand_path(File.join(File.dirname(__FILE__), 'db.bin')), 'w') {|f| f.write(Marshal.dump(@conv_items)) }
    end

    def find_item(inputs, output, input_format=nil, output_format=nil, package_manager=nil)
      inputs = inputs.map{ |input| Pathname(input) }
      output = Pathname(output)
      inputs_ext = inputs.map{ |input| input.extname.slice(1..-1) }
      output_ext = output.extname.slice(1..-1)
      inputs_ext = input_format if input_format
      output_ext = output_format if output_format
      return @conv_items[inputs_ext].select{ |t| t.output == [output_ext] }.sort_by{ |item| item.last_access }.map{ |t| t.input_path=inputs;t.output_path=output;t }
      #@prompt.select
    end

    def list(inputs, input_format=nil, package_manager=nil)
      inputs_path = inputs.map{ |input| Pathname(input) }
      inputs_ext = inputs_path.map{ |input| input.extname.slice(1..-1) }
      inputs_ext = input_format if input_format
      #input = Pathname(input)
      #input_ext = input.extname
      return @conv_items[inputs_ext].sort_by{ |item| item.last_access }.sort_by{ |item| item.last_access }.map{ |t| t.input_path=inputs;t.output_path=t.generate_output_filename(inputs[0]);t }
    end

    private
    def self.rehash
      items = []
      items << ConvItem.new(input: ["bam"], output: ["sam"], command: "samtools view -h {input} > {output}", substitute: "{}", package: "samtools", output_buffer: :stdout)
      items << ConvItem.new(input: ["sam"], output: ["bam"], command: "samtools view -bS {input} > {output}", substitute: "{}", package: "samtools", output_buffer: :stdout)
      items << ConvItem.new(input: ["bam"], output: ["bai"], command: "samtools index {input}", substitute: "{}", package: "samtools")
      items << ConvItem.new(input: ["fa"], output: ["fai"], command: "samtools faidx {input}", substitute: "{}", package: "samtools")
      items << ConvItem.new(input: ["bam"], output: ["fasta"], command: "samtools fasta {input} > {output}", substitute: "{}", package: "samtools")
      items << ConvItem.new(input: ["bam"], output: ["fastq"], command: "samtools fastq {input} > {output}", substitute: "{}", package: "samtools")
      items << ConvItem.new(input: ["sam"], output: ["gtf"], command: "cufflinks {input}", substitute: "{}", package: "cufflinks")
      items << ConvItem.new(input: ["bam"], output: ["gtf"], command: "cufflinks {input}", substitute: "{}", package: "cufflinks")
      items << ConvItem.new(input: ["vcf"], output: ["vcf.gz"], command: "bgzip -c {input} > {output}", substitute: "{}", package: "samtools")
      items << ConvItem.new(input: ["mid"], output: ["wav"], command: "timidity {input} -Ow -o {output}", substitute: "{}", package: "tmidity")
      items << ConvItem.new(input: ["sam", "gtf"], output: ["txt"], command: "htseq-count {input} {input} > {output}", substitute: "{}", package: "htseq", pacman: "pip install")
      items << ConvItem.new(input: ["bam", "gtf"], output: ["txt"], command: "htseq-count -f bam {input} {input} > {output}", substitute: "{}", package: "htseq", pacman: "pip install")
      ["png", "pdf", "ps", "jpg"].each do |i|
        items << ConvItem.new(input: ["dot"], output: [i], command: "dot -T#{i} {input} > {output}", substitute: "{}", package: "graphviz")
      end
      ["png", "jpg", "tif", "bmp", "gif"].permutation(2).each do |i, j|
        items << ConvItem.new(input: [i], output: [j], command: "convert {input} {output}", substitute: "{}", package: "imagemagick")
      end
      items << ConvItem.new(input: ["wav"], output: ["mp3"], command: "lame {input} {output}", substitute: "{}", package: "lame")
      items << ConvItem.new(input: ["wav"], output: ["mp3"], command: "avconv -i {input} {output}", substitute: "{}", package: "libav")
      items << ConvItem.new(input: ["md"], output: ["html"], command: "python -m markdown {input} > {output}", substitute: "{}", package: "markdown", pacman: "pip install")
      items << ConvItem.new(input: ["md"], output: ["html"], command: "pandoc -f markdown_github {input} > {output}", substitute: "{}", package: "pandoc")
      items
    end
  end

  class ConvItem
    attr_accessor :input, :output, :command, :package, :substitute, :last_access, :pacman, :output_buffer, :input_path, :output_path

    def initialize **params
      Config.load_and_set_settings(Dir.home + ".conv")
      @package_manager = Settings.package_manager
      params.each{|k,v| self.send("#{k}=", v) if self.methods.include?(k)}
    end

    def generate_command(options, input, output)
      return [self.cmd, self.brew]
    end

    def generate_output_filename(input, option=:exchange)
      extensions = self.output
      case option
      when :exchange
        extensions = extensions.map{ |ext| Pathname(input).sub_ext("."+ext).to_s }
      when :prefix
        extensions = extensions.map{ |ext| input + "." + ext }
      #when :withtime
      #  extensions = extensions.map{ |ext| Pathname(input + "." + ext) }
      end
      extensions[0]
    end

    def cmd(inputs, output)
      cmd = self.command

      inputs.each{ |input| cmd.sub!("{input}", input)}
      cmd.sub!("{output}", output)
      return cmd
    end

    def to_s
      return self.output.to_s + " " + self.cmd(@input_path, @output_path) + "  (" + self.package + ")"
    end

    def run!(inputs, output)
      raise "Error -- no input file" unless inputs.all?{ |t| Pathname(t).exist?}
      status, stdout, stderr = systemu self.cmd(inputs, output)

      return status, stdout, stderr
    end

    def brew!
      status, stdout, stderr= systemu self.brew_command
      #status, stdout = Open3.capture2(self.brew_command)
      return status, stdout, stderr
    end

    def brew_command
      if self.pacman
        return pacman + " " + package
      else
        if !@package_manager
          return "brew install " + package
        else
          return @package_manager + package
        end
      end
    end
  end

end
