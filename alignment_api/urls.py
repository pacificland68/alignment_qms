from django.urls import path, include
from .views import (
    AlignmentApiView,
    AlignmentNumericalApiView,
    Test
)

urlpatterns = [
    path('api', AlignmentApiView.as_view()),
    path('api/number', AlignmentNumericalApiView.as_view()),
    path('test', Test.as_view())
]