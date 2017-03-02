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

    def find_item(inputs, output, input_format=nil, output_format=nil)
      inputs = inputs.map{ |input| Pathname(input) }
      output = Pathname(output)
      raise "Error -- no input file" unless inputs.all?{ |t| t.exist?}
      inputs_ext = inputs.map{ |input| input.extname.slice(1..-1) }
      output_ext = output.extname.slice(1..-1)
      inputs_ext = input_format if input_format
      output_ext = output_format if output_format
      return @conv_items[inputs_ext].select{ |t| t.output == [output_ext] }
      #@prompt.select
    end

    def list(inputs, input_format=nil)
      inputs = inputs.map{ |input| Pathname(input) }
      inputs_ext = inputs.map{ |input| input.extname.slice(1..-1) }
      inputs_ext = input_format if input_format
      #input = Pathname(input)
      #input_ext = input.extname
      return @conv_items[inputs_ext]
    end

    private
    def rehash
      items = []
      items << ConvItem.new(input: ["bam"], output: ["sam"], command: "samtools view -bS {input} > {output}", substitute: "{}", package: "samtools", output_buffer: :stdout)
      items << ConvItem.new(input: ["sam"], output: ["bam"], command: "samtools view -h {input} > {output}", substitute: "{}", package: "samtools", output_buffer: :stdout)
      items << ConvItem.new(input: ["bam"], output: ["bai"], command: "samtools index {input}", substitute: "{}", package: "samtools")
      items << ConvItem.new(input: ["fa"], output: ["fai"], command: "samtools faidx {input}", substitute: "{}", package: "samtools")
      items << ConvItem.new(input: ["dot"], output: ["png"], command: "dot -Tpng {input} > {output}", substitute: "{}", package: "graphviz")
      items << ConvItem.new(input: ["png"], output: ["jpg"], command: "convert -format jpg {input} > {output}", substitute: "{}", package: "imagemagick")
      items << ConvItem.new(input: ["md"], output: ["html"], command: "python -m markdown {input} > {output}", substitute: "{}", package: "markdown", pacman: "pip install")
      items
    end
  end

  class ConvItem
    attr_accessor :input, :output, :command, :package, :substitute, :priority, :pacman, :output_buffer

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
      return self.output.to_s + " " + self.command
    end

    def run!(inputs, output)
      status, stdout, stderr = systemu self.cmd(inputs, output)

      return status, stdout, stderr
    end

    def brew!
      status, stdout, stderr = systemu self.brew_command
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
