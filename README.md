# FFC
[![Gem Version](https://badge.fury.io/rb/ffc.svg)](https://badge.fury.io/rb/ffc)
[![Build Status](https://travis-ci.org/6br/ffc.svg?branch=master)](https://travis-ci.org/6br/ffc)
[![Stories in Ready](https://badge.waffle.io/6br/ffc.svg?label=ready&title=Ready)](http://waffle.io/6br/ffc)
[![Code Climate](https://codeclimate.com/github/6br/ffc/badges/gpa.svg)](https://codeclimate.com/github/6br/ffc)

File Formats Conversion Tool

Especially, this tool may be useful for bioinformatics.

## Installation

Install it yourself as:

    $ gem install ffc

## Examples

To compile programs:

    $ ffc convert sample.c a.out
    $ ffc c sample.c a.out

To convert markdown file to HTML file:

    $ ffc convert a.md a.html
    $ ffc c a.md a.html

To convert "sam" format file to "bam" format file (It converts automatically if you do not know appreciate commands and parameters.):

    $ ffc convert a.sam a.bam
    $ ffc c a.sam a.bam
    
To show the whole list of available convert commands:

    $ ffc list a.png
    $ ffc l a.png
    
To show help:

    $ ffc help [COMMAND]
    $ ffc h [COMMAND]

When there are alternative tools, you can select one of them.

## Development

After checking out the repo, run `bin/setup` to install dependencies. Then, run `rake spec` to run the tests. You can also run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release`, which will create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).

## Notice

If there are no package to run a suggested command, this tool will install the package with `brew install` in default.

If you use another package manager such as `apt`, `yum` or `conda`, you can configure the default command. But there is no interface or guide yet. 

## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/6br/ffc.

