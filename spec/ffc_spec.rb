require "spec_helper"

RSpec.describe Conv do
  it "has a version number" do
    expect(Conv::VERSION).not_to be nil
  end

  it "should display output" do
    expect {Conv::CLI.new.invoke(:convert, ["spec/sample/a.sam", "b.non"], {hidelist: true})}.to output("\e[31mNo candidates.\e[0m\n").to_stdout
  end

end

