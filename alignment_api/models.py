from django.db import models
from django.contrib.auth.models import User

class alignment(models.Model):
    task = models.CharField(max_length=180)

    def __self__(self):
        return self.task

# Create your models here.
