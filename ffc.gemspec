# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'ffc/version'

Gem::Specification.new do |spec|
  spec.name          = "ffc"
  spec.version       = Conv::VERSION
  spec.authors       = ["6br"]
  spec.email         = ["git_hub-ruby@yahoo.co.jp"]

  spec.summary       = %q{Support to convert file formats by preset UNIX commands}
  spec.description   = %q{Support to convert file formats choosing from a preset UNIX command list}
  spec.homepage      = "https://github.com/6br/ffc"

  # Prevent pushing this gem to RubyGems.org. To allow pushes either set the 'allowed_push_host'
  # to allow pushing to a single host or delete this section to allow pushing to any host.
  if spec.respond_to?(:metadata)
    #spec.metadata['allowed_push_host'] = "http://mygemserver.com"
  else
    raise "RubyGems 2.0 or newer is required to protect against " \
      "public gem pushes."
  end

  spec.files         = `git ls-files -z`.split("\x0").reject do |f|
    f.match(%r{^(test|spec|features)/})
  end
  spec.bindir        = "exe"
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler"
  spec.add_development_dependency "rake", ">= 12.3.3"
  spec.add_development_dependency "rspec", "~> 3.0"
  spec.add_development_dependency "coveralls"
  spec.add_development_dependency "codeclimate-test-reporter"
  spec.add_runtime_dependency "thor"
  spec.add_runtime_dependency "tty-prompt"
  spec.add_runtime_dependency "config"
  spec.add_runtime_dependency "systemu"
end
