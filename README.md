# FFC

File Formats Convertion Tool

Especially, this tool may be useful for bioinformatics.

## Installation

Install it yourself as:

    $ gem install ffc

## Examples

* To convert markdown file to html file

    $ ffc convert a.md a.html
    $ ffc c a.md a.html
    
* To show the whole list of inputs format

    $ ffc list a.png
    $ ffc l a.png

When there are alternative tools, you can select one of them.

## Development

After checking out the repo, run `bin/setup` to install dependencies. Then, run `rake spec` to run the tests. You can also run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release`, which will create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).

## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/6br/ffc.

