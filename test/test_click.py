import click

@click.command()
@click.option('-w', '--width', type=int, default=0)
@click.option('--option2')
@click.argument('argument')
def app(width, option2, argument):
    click.echo("params: {} {} {}".format(width, option2, argument))


def test_app():

    print ('now')

    # app(["-w" 3, ])
    print ('run')
    # print (app(["arg", "--option2", "4", "-w", 3]))

    # app(["arg", "-w", 3, "--option2", "4"])
    #
    app(["-w", 3, "--option2", "4", "arg"])
