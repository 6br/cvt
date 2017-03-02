# -*- coding: utf-8 -*-
require 'ffc'
require 'thor'
require 'tty-prompt'

module Conv
  class CLI < Thor
    package_name "ffc"
    #default_command :convert

    desc "convert INPUTS OUTPUT", "convert file formats from input to output."
    option :confirm, type: :boolean, aliases: '-c', desc: "Confirm or choice from alternative tools."
    #option :onetime, type: :boolean, aliases: '-1', desc: "Do not save in the internal command history."
    option :overwrite, type: :boolean, aliases: '-w', desc: "Overwrite output file without confirm."
    option :autobrew, type: :boolean, aliases: '-b', desc: "Install without confirmation if there is no candidate binary."
    option :input_format, type: :array, aliases: '-i', desc: "Specify input file formats regardless of extensions."
    option :output_format, type: :string, aliases: '-o', desc: "Specify the output file format regardless of the extension."
    option :hidelist, type: :boolean, aliases: '-h', desc: "Do not show lists related to inputs when no candidates"
    option :quiet, type: :boolean, aliases: '-q', desc: "Run quietly."
    option :dryrun, type: :boolean, aliases: '-d', desc: "Do not execute command."
    option :package_manager, type: :array, aliases: '-p', desc: "Set Package manager prefix such as 'brew install'"
    def convert(*inputs, output) # Define commands as a method.
      prompt = TTY::Prompt.new(interrupt: :exit)
      items = Conv.new.find_item(inputs, output, input_format=input_format, output_format=output_format, package_manager: options[:package_manager])
      candidate = items[0]
      if items.empty?
        prompt.error "No candidates."
        inputs << output
        list(*inputs) if !options[:hidelist]
      else
          if options[:confirm]#interactive]
            candidate = prompt.select("Choose commands (To cancel, press CTRL+C)", items) #do |item|
          end
          prompt.ok candidate.cmd(inputs, output) unless options[:quiet]
          return if !options[:overwrite] && !ask_overwrite(output) || options[:dryrun]
          status, stdout, stderr = candidate.run!(inputs, output)
          prompt.ok stdout
          prompt.error stderr
          if status != 0
            if options[:autobrew] || !prompt.no?("Install " + candidate.package + " with package manager? (" + candidate.brew_command + ")")
              status, stdout, stderr = candidate.brew!
              prompt.say stdout
              prompt.error stderr
              if status == 0
                status, stdout, stderr = candidate.run!(inputs, output)
                prompt.ok stdout
                prompt.error stderr
              end
            end
          else
            Conv.update!(candidate) unless options[:onetime]
          end
      end
    end

    desc "list INPUTS", "show aveilable convert formats."
    #option :interactive, type: :boolean, aliases: '-i', desc: "Only show lists"
    option :showonly, type: :boolean, aliases: '-s', desc: "Only show lists"
    option :add_prefix, type: :boolean, aliases: '-p', desc: "Add prefix to input's filename"
    option :overwrite, type: :boolean, aliases: '-w', desc: "Overwrite output without confirm"
    option :autobrew, type: :boolean, aliases: '-b', desc: "Install without confirmation if there is no candidate binary."
    option :input_format, type: :array, aliases: '-i', desc: "Specify input file formats regardless of extensions."
    option :quiet, type: :boolean, aliases: '-q', desc: "Run quietly."
    option :dryrun, type: :boolean, aliases: '-d', desc: "Do not execute command."
    option :package_manager, type: :array, aliases: '-p', desc: "Set Package manager prefix such as 'brew install'"
    def list(*inputs)
      prompt = TTY::Prompt.new(interrupt: :exit)
      items = Conv.new.list(inputs, input_format=input_format, package_manager: options[:package_manager])
      #p inputs,items
      if items.empty?
        prompt.error "No candidates."
        return
      end
      if !options[:showonly]
        candidate = prompt.select("Choose commands (To cancel: CTRL+C)", items) #do |item|
        option = options[:add_prefix] ? :prefix : :exchange
        output = candidate.generate_output_filename(inputs[0], option)

        return if !options[:overwrite] && !ask_overwrite(output) || options[:dryrun]
        prompt.ok candidate.cmd(inputs, output) unless options[:quiet]
        status, stdout, stderr = candidate.run!(inputs, output)#(options[:autobrew])
        prompt.say stdout
        prompt.error stderr
        if status != 0
          if options[:autobrew] || !prompt.no?("Install " + candidate.package + " with package manager? (" + candidate.brew_command + ")")
              status, stdout, stderr = candidate.brew!
              prompt.say stdout
              prompt.error stderr
              if status == 0
                status, stdout, stderr = candidate.run!(inputs, output)
                prompt.say stdout
                prompt.error stderr
              end
          end
        else
          Conv.update!(candidate) unless options[:onetime]
        end
      else
        prompt.say items.join("")

      end
    end

    desc "sort INPUT", "Sort the file. (Not Implemented)"
    option :overwrite, type: :boolean, aliases: '-w', desc: "Overwrite existing file"
    def sort(input)
#
    end

    desc "index INPUT", "Index the file. (Not Implemented)"
    option :overwrite, type: :boolean, aliases: '-w', desc: "Overwrite existing file"
    def index(input)
#
    end

    desc "rehash", "Rehash index and clean command history. (For debugging)"
    def rehash
      Conv.save!
    end

    private
    def ask_overwrite(output)
      if Pathname(output).exist?
        return TTY::Prompt.new.yes?("Overwrite existing file " + output + "?")
      end
      true
    end

  end
end

