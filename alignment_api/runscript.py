from django.core.management.base import BaseCommand

class Command(BaseCommand):
    help = 'Runs a custom script in the Django environment.'

    def handle(self, *args, **options):
        # 下面是你的 Python 脚本代码
        print("Hello from the custom script!")
        # ... 更多代码 ...
