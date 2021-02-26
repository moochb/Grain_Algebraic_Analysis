import scipy.special
from itertools import combinations 
import time
import random

startprecomp = time.perf_counter()
#define the paramaters

#how many clocks are required
no_clocks = 15000

#the size of the register
state_size = 20

#how many output bits are needed
no_output_bits=15000

#degree of the system of equations to be solved
degree = 2

#the number of monomials of degree 1->max degree
#the variables 'no_monomials' will be the number of linearisation variables
#'scipy.special.comb' provides the binomial coefficient

s_out = [0]*state_size
b_out = [0]*state_size

for i in range(state_size):
    s_out[i] = random.randint(0,1)
    b_out[i] = random.randint(0,1)
    
print(s_out)
print(b_out)

add=[0]*no_clocks
S_out = s_out+add
B_out = b_out+add



O = [0]*(no_clocks-state_size)
for i in range(state_size,no_clocks):
    O[i-state_size] =mod(S_out[i-(state_size-4)]*S_out[i-(state_size-1)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-4)]*S_out[i-(state_size-9)]*B_out[i-(state_size-12 )]+ S_out[i-(state_size-4 )]+ S_out[i-(state_size-1)]*S_out[i-(state_size-13)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-1)]*S_out[i-(state_size-9)]*B_out[i-(state_size-12 )]+ S_out[i-(state_size-1 )]+ S_out[i-(state_size-13)]*S_out[i-(state_size-9)]*B_out[i-(state_size-12 )]+ S_out[i-(state_size-13)]*S_out[i-(state_size-9 )]+ S_out[i-(state_size-13)]*B_out[i-(state_size-12 )]+ S_out[i-(state_size-13 )]+ B_out[i-(state_size-12)],2)
    B_out[i] = mod(S_out[i-(state_size-0 )]+ B_out[i-(state_size-0 )]+ B_out[i-(state_size-13 )]+ B_out[i-(state_size-19 )]+ B_out[i-(state_size-15 )]+ B_out[i-(state_size-2)]*B_out[i-(state_size-15 )]+ B_out[i-(state_size-3)]*B_out[i-(state_size-5 )]+ B_out[i-(state_size-7)]*B_out[i-(state_size-8 )]+ B_out[i-(state_size-14)]*B_out[i-(state_size-19 )]+ B_out[i-(state_size-13)]*B_out[i-(state_size-6)]*B_out[i-(state_size-17)]*B_out[i-(state_size-18 )]+ B_out[i-(state_size-10)]*B_out[i-(state_size-11)]*B_out[i-(state_size-12 )],2)
    S_out[i] = mod((S_out[i-20]+S_out[i-9]+S_out[i-5]+S_out[i-3]),2)


no_monomials = 0
for i in range(1,degree+1):
    no_monomials = no_monomials + binomial(state_size,i)

#define polynomial ring in which the the INITIAL STATE, OUTPUT and LINEARISATION variables will be defined
R=BooleanPolynomialRing(state_size+no_monomials,'x')
#Inject these variables for use
R.inject_variables()

#define initial state
s = [0]*state_size
for i in range(state_size):
    s[i] = eval("x" + str(i))

# print(s)
#define output
# Y = [0]*no_output_bits
# for i in range(no_output_bits):
#     Y[i] = eval("x" + str(state_size+i))

# print(Y)
#We can consider 's' as a register that will be shifted for a specified number of times. Thus, there will exist more entries than simply those defined in 's'.
#We therefore initialise a vector 'add' that will be appended to s to hold all the update bits (one vector for equations with respect to keystream bits S and one with respect the the linear approximation$
Alpha = [6195 ,  6194 ,  6193 ,  6191 ,  6190 ,  6189 ,  6187 ,  6186 ,  6185 ,  6184 ,  6176 ,  6174 ,  6173 ,  6172 ,  6169 ,  6167 ,  6166 ,  6163 ,  6162 ,  6161 ,  6160 ,  6159 ,  6158 ,  6157 ,  6156 ,  6155 ,  6151 ,  6150 ,  6147 ,  6144 ,  6141 ,  6139 ,  6138 ,  6137 ,  6136 ,  6135 ,  6134 ,  6133 ,  6132 ,  6130 ,  6128 ,  6127 ,  6125 ,  6123 ,  6121 ,  6120 ,  6117 ,  6115 ,  6114 ,  6111 ,  6108 ,  6107 ,  6104 ,  6103 ,  6099 ,  6096 ,  6095 ,  6093 ,  6092 ,  6091 ,  6087 ,  6086 ,  6082 ,  6080 ,  6075 ,  6071 ,  6070 ,  6069 ,  6063 ,  6062 ,  6061 ,  6058 ,  6056 ,  6053 ,  6052 ,  6050 ,  6049 ,  6046 ,  6043 ,  6041 ,  6038 ,  6036 ,  6035 ,  6033 ,  6029 ,  6028 ,  6027 ,  6025 ,  6024 ,  6018 ,  6016 ,  6014 ,  6010 ,  6005 ,  6004 ,  6002 ,  6001 ,  6000 ,  5998 ,  5997 ,  5996 ,  5995 ,  5993 ,  5991 ,  5990 ,  5989 ,  5986 ,  5985 ,  5984 ,  5978 ,  5977 ,  5975 ,  5972 ,  5970 ,  5967 ,  5963 ,  5960 ,  5959 ,  5956 ,  5955 ,  5954 ,  5947 ,  5946 ,  5945 ,  5943 ,  5941 ,  5940 ,  5938 ,  5937 ,  5934 ,  5933 ,  5928 ,  5924 ,  5922 ,  5921 ,  5912 ,  5911 ,  5910 ,  5909 ,  5908 ,  5906 ,  5905 ,  5904 ,  5901 ,  5900 ,  5897 ,  5896 ,  5895 ,  5894 ,  5893 ,  5890 ,  5888 ,  5887 ,  5885 ,  5882 ,  5881 ,  5880 ,  5878 ,  5876 ,  5875 ,  5874 ,  5873 ,  5872 ,  5870 ,  5869 ,  5866 ,  5865 ,  5861 ,  5858 ,  5855 ,  5853 ,  5850 ,  5849 ,  5843 ,  5842 ,  5840 ,  5838 ,  5835 ,  5831 ,  5830 ,  5828 ,  5826 ,  5825 ,  5824 ,  5823 ,  5822 ,  5819 ,  5818 ,  5817 ,  5814 ,  5812 ,  5811 ,  5808 ,  5806 ,  5805 ,  5803 ,  5800 ,  5798 ,  5797 ,  5795 ,  5794 ,  5793 ,  5792 ,  5790 ,  5789 ,  5787 ,  5782 ,  5780 ,  5776 ,  5775 ,  5774 ,  5772 ,  5770 ,  5768 ,  5767 ,  5764 ,  5761 ,  5760 ,  5759 ,  5758 ,  5756 ,  5753 ,  5749 ,  5748 ,  5747 ,  5744 ,  5742 ,  5740 ,  5738 ,  5736 ,  5731 ,  5730 ,  5729 ,  5727 ,  5726 ,  5725 ,  5721 ,  5718 ,  5716 ,  5714 ,  5711 ,  5710 ,  5709 ,  5705 ,  5702 ,  5699 ,  5697 ,  5696 ,  5695 ,  5693 ,  5691 ,  5687 ,  5686 ,  5684 ,  5683 ,  5682 ,  5680 ,  5679 ,  5676 ,  5675 ,  5674 ,  5672 ,  5671 ,  5670 ,  5669 ,  5668 ,  5666 ,  5660 ,  5659 ,  5658 ,  5657 ,  5650 ,  5648 ,  5645 ,  5643 ,  5638 ,  5637 ,  5631 ,  5630 ,  5629 ,  5627 ,  5625 ,  5624 ,  5621 ,  5619 ,  5617 ,  5614 ,  5613 ,  5612 ,  5611 ,  5609 ,  5607 ,  5606 ,  5605 ,  5604 ,  5602 ,  5601 ,  5600 ,  5598 ,  5595 ,  5593 ,  5592 ,  5591 ,  5587 ,  5586 ,  5584 ,  5583 ,  5581 ,  5580 ,  5575 ,  5574 ,  5573 ,  5571 ,  5570 ,  5569 ,  5568 ,  5566 ,  5563 ,  5562 ,  5560 ,  5559 ,  5556 ,  5553 ,  5552 ,  5548 ,  5547 ,  5546 ,  5542 ,  5540 ,  5539 ,  5538 ,  5537 ,  5533 ,  5532 ,  5529 ,  5527 ,  5523 ,  5522 ,  5521 ,  5520 ,  5515 ,  5514 ,  5513 ,  5512 ,  5507 ,  5506 ,  5504 ,  5501 ,  5498 ,  5495 ,  5493 ,  5491 ,  5488 ,  5487 ,  5483 ,  5481 ,  5480 ,  5473 ,  5472 ,  5466 ,  5464 ,  5463 ,  5461 ,  5460 ,  5454 ,  5450 ,  5449 ,  5448 ,  5446 ,  5444 ,  5442 ,  5440 ,  5438 ,  5435 ,  5433 ,  5432 ,  5431 ,  5430 ,  5429 ,  5428 ,  5424 ,  5419 ,  5418 ,  5416 ,  5414 ,  5413 ,  5412 ,  5410 ,  5407 ,  5406 ,  5404 ,  5403 ,  5402 ,  5401 ,  5400 ,  5396 ,  5394 ,  5390 ,  5389 ,  5388 ,  5387 ,  5385 ,  5384 ,  5382 ,  5381 ,  5378 ,  5377 ,  5376 ,  5375 ,  5374 ,  5371 ,  5368 ,  5367 ,  5366 ,  5365 ,  5363 ,  5360 ,  5357 ,  5355 ,  5353 ,  5351 ,  5350 ,  5348 ,  5347 ,  5346 ,  5345 ,  5344 ,  5343 ,  5338 ,  5335 ,  5332 ,  5331 ,  5330 ,  5328 ,  5318 ,  5317 ,  5316 ,  5315 ,  5313 ,  5311 ,  5306 ,  5305 ,  5302 ,  5301 ,  5300 ,  5299 ,  5296 ,  5295 ,  5292 ,  5290 ,  5288 ,  5286 ,  5285 ,  5283 ,  5281 ,  5278 ,  5273 ,  5272 ,  5271 ,  5266 ,  5264 ,  5262 ,  5261 ,  5259 ,  5258 ,  5257 ,  5254 ,  5253 ,  5252 ,  5250 ,  5249 ,  5247 ,  5246 ,  5242 ,  5241 ,  5240 ,  5239 ,  5238 ,  5237 ,  5236 ,  5233 ,  5232 ,  5229 ,  5228 ,  5223 ,  5219 ,  5216 ,  5215 ,  5214 ,  5213 ,  5210 ,  5209 ,  5207 ,  5203 ,  5202 ,  5201 ,  5200 ,  5198 ,  5197 ,  5195 ,  5194 ,  5193 ,  5192 ,  5190 ,  5189 ,  5187 ,  5186 ,  5183 ,  5181 ,  5180 ,  5179 ,  5178 ,  5176 ,  5174 ,  5173 ,  5171 ,  5170 ,  5169 ,  5166 ,  5164 ,  5162 ,  5155 ,  5153 ,  5152 ,  5151 ,  5149 ,  5147 ,  5146 ,  5144 ,  5142 ,  5141 ,  5140 ,  5135 ,  5134 ,  5133 ,  5131 ,  5129 ,  5126 ,  5125 ,  5123 ,  5122 ,  5120 ,  5118 ,  5116 ,  5115 ,  5113 ,  5112 ,  5111 ,  5109 ,  5106 ,  5104 ,  5102 ,  5101 ,  5099 ,  5096 ,  5094 ,  5091 ,  5090 ,  5087 ,  5085 ,  5084 ,  5081 ,  5080 ,  5078 ,  5076 ,  5075 ,  5074 ,  5070 ,  5068 ,  5066 ,  5065 ,  5063 ,  5062 ,  5060 ,  5058 ,  5053 ,  5049 ,  5048 ,  5046 ,  5042 ,  5039 ,  5038 ,  5036 ,  5033 ,  5030 ,  5028 ,  5027 ,  5026 ,  5024 ,  5023 ,  5022 ,  5020 ,  5019 ,  5015 ,  5011 ,  5007 ,  5006 ,  5005 ,  5003 ,  5002 ,  5001 ,  5000 ,  4996 ,  4992 ,  4991 ,  4987 ,  4984 ,  4982 ,  4980 ,  4979 ,  4978 ,  4977 ,  4976 ,  4975 ,  4971 ,  4968 ,  4965 ,  4960 ,  4958 ,  4957 ,  4955 ,  4954 ,  4953 ,  4952 ,  4950 ,  4944 ,  4943 ,  4941 ,  4940 ,  4939 ,  4938 ,  4936 ,  4935 ,  4934 ,  4930 ,  4929 ,  4928 ,  4927 ,  4924 ,  4923 ,  4921 ,  4918 ,  4917 ,  4914 ,  4911 ,  4907 ,  4906 ,  4904 ,  4903 ,  4900 ,  4898 ,  4897 ,  4895 ,  4894 ,  4893 ,  4891 ,  4889 ,  4888 ,  4886 ,  4883 ,  4882 ,  4881 ,  4877 ,  4873 ,  4872 ,  4871 ,  4867 ,  4865 ,  4864 ,  4862 ,  4861 ,  4859 ,  4857 ,  4856 ,  4852 ,  4851 ,  4849 ,  4848 ,  4847 ,  4843 ,  4839 ,  4838 ,  4835 ,  4834 ,  4833 ,  4832 ,  4829 ,  4827 ,  4826 ,  4824 ,  4823 ,  4822 ,  4814 ,  4810 ,  4806 ,  4804 ,  4800 ,  4793 ,  4791 ,  4790 ,  4787 ,  4786 ,  4783 ,  4782 ,  4779 ,  4777 ,  4775 ,  4773 ,  4772 ,  4770 ,  4769 ,  4768 ,  4766 ,  4764 ,  4760 ,  4756 ,  4754 ,  4753 ,  4751 ,  4750 ,  4748 ,  4745 ,  4742 ,  4737 ,  4735 ,  4733 ,  4732 ,  4731 ,  4728 ,  4726 ,  4724 ,  4723 ,  4718 ,  4717 ,  4712 ,  4710 ,  4709 ,  4707 ,  4702 ,  4700 ,  4699 ,  4698 ,  4696 ,  4694 ,  4692 ,  4688 ,  4686 ,  4685 ,  4681 ,  4679 ,  4678 ,  4673 ,  4670 ,  4669 ,  4665 ,  4664 ,  4662 ,  4659 ,  4653 ,  4650 ,  4649 ,  4648 ,  4647 ,  4643 ,  4642 ,  4639 ,  4637 ,  4635 ,  4634 ,  4633 ,  4630 ,  4629 ,  4628 ,  4627 ,  4623 ,  4622 ,  4621 ,  4613 ,  4612 ,  4609 ,  4606 ,  4604 ,  4603 ,  4601 ,  4594 ,  4592 ,  4588 ,  4586 ,  4584 ,  4582 ,  4578 ,  4577 ,  4575 ,  4574 ,  4571 ,  4570 ,  4569 ,  4568 ,  4566 ,  4564 ,  4563 ,  4562 ,  4561 ,  4558 ,  4557 ,  4555 ,  4554 ,  4553 ,  4552 ,  4550 ,  4549 ,  4546 ,  4542 ,  4541 ,  4539 ,  4538 ,  4537 ,  4536 ,  4533 ,  4530 ,  4529 ,  4525 ,  4521 ,  4517 ,  4516 ,  4514 ,  4512 ,  4510 ,  4507 ,  4506 ,  4505 ,  4501 ,  4500 ,  4499 ,  4498 ,  4497 ,  4496 ,  4489 ,  4488 ,  4486 ,  4484 ,  4483 ,  4481 ,  4479 ,  4478 ,  4477 ,  4476 ,  4474 ,  4465 ,  4464 ,  4460 ,  4457 ,  4456 ,  4448 ,  4446 ,  4445 ,  4441 ,  4437 ,  4435 ,  4434 ,  4433 ,  4432 ,  4431 ,  4429 ,  4427 ,  4425 ,  4423 ,  4419 ,  4418 ,  4417 ,  4415 ,  4414 ,  4413 ,  4412 ,  4409 ,  4408 ,  4406 ,  4402 ,  4401 ,  4400 ,  4398 ,  4397 ,  4393 ,  4392 ,  4391 ,  4390 ,  4389 ,  4387 ,  4386 ,  4383 ,  4381 ,  4379 ,  4377 ,  4374 ,  4372 ,  4370 ,  4369 ,  4368 ,  4363 ,  4362 ,  4361 ,  4360 ,  4356 ,  4354 ,  4353 ,  4352 ,  4351 ,  4350 ,  4349 ,  4348 ,  4346 ,  4344 ,  4343 ,  4342 ,  4339 ,  4337 ,  4335 ,  4334 ,  4333 ,  4332 ,  4331 ,  4330 ,  4329 ,  4327 ,  4326 ,  4325 ,  4324 ,  4323 ,  4321 ,  4320 ,  4315 ,  4312 ,  4310 ,  4307 ,  4306 ,  4304 ,  4302 ,  4300 ,  4299 ,  4298 ,  4295 ,  4294 ,  4292 ,  4287 ,  4283 ,  4282 ,  4279 ,  4278 ,  4277 ,  4276 ,  4274 ,  4270 ,  4266 ,  4262 ,  4260 ,  4259 ,  4258 ,  4256 ,  4255 ,  4254 ,  4249 ,  4247 ,  4246 ,  4240 ,  4239 ,  4238 ,  4237 ,  4235 ,  4232 ,  4228 ,  4226 ,  4225 ,  4224 ,  4222 ,  4221 ,  4220 ,  4218 ,  4216 ,  4215 ,  4209 ,  4206 ,  4203 ,  4202 ,  4198 ,  4197 ,  4195 ,  4194 ,  4193 ,  4192 ,  4191 ,  4187 ,  4186 ,  4184 ,  4183 ,  4181 ,  4179 ,  4176 ,  4175 ,  4169 ,  4168 ,  4167 ,  4164 ,  4162 ,  4161 ,  4160 ,  4159 ,  4158 ,  4156 ,  4152 ,  4151 ,  4150 ,  4149 ,  4148 ,  4147 ,  4146 ,  4145 ,  4144 ,  4140 ,  4139 ,  4136 ,  4135 ,  4133 ,  4132 ,  4130 ,  4129 ,  4127 ,  4126 ,  4124 ,  4123 ,  4122 ,  4121 ,  4115 ,  4114 ,  4113 ,  4112 ,  4111 ,  4110 ,  4109 ,  4108 ,  4106 ,  4105 ,  4104 ,  4102 ,  4101 ,  4100 ,  4098 ,  4097 ,  4096 ,  4095 ,  4094 ,  4093 ,  4086 ,  4084 ,  4083 ,  4082 ,  4081 ,  4078 ,  4076 ,  4074 ,  4073 ,  4070 ,  4068 ,  4067 ,  4066 ,  4064 ,  4063 ,  4062 ,  4061 ,  4060 ,  4057 ,  4056 ,  4054 ,  4052 ,  4048 ,  4047 ,  4046 ,  4044 ,  4039 ,  4033 ,  4031 ,  4029 ,  4028 ,  4025 ,  4022 ,  4020 ,  4019 ,  4018 ,  4017 ,  4012 ,  4011 ,  4009 ,  4008 ,  4007 ,  4006 ,  4005 ,  4004 ,  4003 ,  4002 ,  4001 ,  4000 ,  3999 ,  3998 ,  3995 ,  3991 ,  3988 ,  3984 ,  3983 ,  3977 ,  3974 ,  3971 ,  3968 ,  3965 ,  3964 ,  3963 ,  3962 ,  3961 ,  3958 ,  3956 ,  3955 ,  3951 ,  3949 ,  3946 ,  3945 ,  3943 ,  3939 ,  3938 ,  3934 ,  3933 ,  3929 ,  3928 ,  3926 ,  3924 ,  3921 ,  3917 ,  3916 ,  3915 ,  3914 ,  3913 ,  3912 ,  3911 ,  3904 ,  3903 ,  3902 ,  3900 ,  3899 ,  3898 ,  3897 ,  3896 ,  3894 ,  3890 ,  3889 ,  3885 ,  3882 ,  3881 ,  3879 ,  3877 ,  3876 ,  3875 ,  3874 ,  3871 ,  3870 ,  3869 ,  3868 ,  3866 ,  3865 ,  3861 ,  3858 ,  3857 ,  3856 ,  3854 ,  3853 ,  3852 ,  3851 ,  3850 ,  3848 ,  3845 ,  3843 ,  3842 ,  3836 ,  3834 ,  3832 ,  3831 ,  3830 ,  3829 ,  3828 ,  3827 ,  3826 ,  3825 ,  3817 ,  3815 ,  3814 ,  3812 ,  3811 ,  3810 ,  3809 ,  3808 ,  3806 ,  3804 ,  3802 ,  3799 ,  3796 ,  3795 ,  3792 ,  3789 ,  3787 ,  3785 ,  3784 ,  3782 ,  3777 ,  3773 ,  3772 ,  3771 ,  3768 ,  3767 ,  3763 ,  3761 ,  3758 ,  3756 ,  3755 ,  3752 ,  3741 ,  3738 ,  3736 ,  3734 ,  3732 ,  3731 ,  3728 ,  3726 ,  3725 ,  3724 ,  3721 ,  3720 ,  3717 ,  3716 ,  3715 ,  3714 ,  3712 ,  3711 ,  3709 ,  3706 ,  3705 ,  3701 ,  3699 ,  3698 ,  3695 ,  3693 ,  3692 ,  3691 ,  3689 ,  3687 ,  3686 ,  3685 ,  3679 ,  3677 ,  3676 ,  3674 ,  3673 ,  3668 ,  3665 ,  3663 ,  3661 ,  3659 ,  3658 ,  3657 ,  3656 ,  3654 ,  3653 ,  3652 ,  3651 ,  3650 ,  3648 ,  3647 ,  3643 ,  3641 ,  3634 ,  3633 ,  3632 ,  3631 ,  3630 ,  3629 ,  3626 ,  3625 ,  3624 ,  3619 ,  3618 ,  3617 ,  3615 ,  3614 ,  3612 ,  3608 ,  3607 ,  3606 ,  3603 ,  3601 ,  3600 ,  3598 ,  3595 ,  3594 ,  3593 ,  3591 ,  3590 ,  3589 ,  3588 ,  3584 ,  3583 ,  3580 ,  3577 ,  3576 ,  3573 ,  3572 ,  3571 ,  3570 ,  3567 ,  3566 ,  3564 ,  3563 ,  3561 ,  3559 ,  3558 ,  3557 ,  3554 ,  3552 ,  3547 ,  3544 ,  3541 ,  3540 ,  3539 ,  3535 ,  3531 ,  3530 ,  3529 ,  3528 ,  3520 ,  3516 ,  3514 ,  3513 ,  3512 ,  3511 ,  3508 ,  3507 ,  3504 ,  3503 ,  3502 ,  3499 ,  3498 ,  3497 ,  3491 ,  3489 ,  3486 ,  3484 ,  3482 ,  3480 ,  3474 ,  3472 ,  3471 ,  3470 ,  3467 ,  3465 ,  3463 ,  3462 ,  3459 ,  3457 ,  3456 ,  3455 ,  3453 ,  3447 ,  3443 ,  3441 ,  3440 ,  3439 ,  3438 ,  3432 ,  3431 ,  3430 ,  3425 ,  3424 ,  3423 ,  3421 ,  3420 ,  3416 ,  3415 ,  3414 ,  3412 ,  3410 ,  3409 ,  3407 ,  3406 ,  3405 ,  3403 ,  3401 ,  3400 ,  3398 ,  3397 ,  3396 ,  3393 ,  3392 ,  3391 ,  3390 ,  3387 ,  3386 ,  3385 ,  3384 ,  3383 ,  3381 ,  3379 ,  3378 ,  3377 ,  3376 ,  3374 ,  3372 ,  3371 ,  3370 ,  3369 ,  3368 ,  3366 ,  3365 ,  3361 ,  3360 ,  3359 ,  3356 ,  3351 ,  3350 ,  3348 ,  3344 ,  3342 ,  3341 ,  3340 ,  3339 ,  3337 ,  3334 ,  3333 ,  3331 ,  3330 ,  3327 ,  3324 ,  3323 ,  3322 ,  3321 ,  3319 ,  3318 ,  3317 ,  3316 ,  3315 ,  3314 ,  3313 ,  3312 ,  3306 ,  3301 ,  3298 ,  3297 ,  3296 ,  3294 ,  3291 ,  3290 ,  3289 ,  3288 ,  3287 ,  3284 ,  3282 ,  3281 ,  3280 ,  3278 ,  3275 ,  3273 ,  3272 ,  3265 ,  3262 ,  3261 ,  3260 ,  3259 ,  3257 ,  3256 ,  3253 ,  3251 ,  3250 ,  3248 ,  3242 ,  3239 ,  3237 ,  3235 ,  3234 ,  3230 ,  3229 ,  3228 ,  3225 ,  3223 ,  3222 ,  3219 ,  3218 ,  3217 ,  3216 ,  3212 ,  3211 ,  3210 ,  3208 ,  3207 ,  3206 ,  3205 ,  3204 ,  3202 ,  3201 ,  3198 ,  3197 ,  3196 ,  3195 ,  3194 ,  3193 ,  3191 ,  3190 ,  3189 ,  3188 ,  3187 ,  3185 ,  3184 ,  3182 ,  3180 ,  3179 ,  3178 ,  3174 ,  3172 ,  3165 ,  3163 ,  3162 ,  3161 ,  3160 ,  3159 ,  3158 ,  3157 ,  3149 ,  3147 ,  3145 ,  3143 ,  3140 ,  3139 ,  3136 ,  3133 ,  3132 ,  3131 ,  3129 ,  3128 ,  3126 ,  3124 ,  3123 ,  3122 ,  3121 ,  3120 ,  3116 ,  3112 ,  3108 ,  3106 ,  3104 ,  3103 ,  3102 ,  3096 ,  3095 ,  3094 ,  3093 ,  3091 ,  3087 ,  3086 ,  3085 ,  3082 ,  3081 ,  3080 ,  3078 ,  3077 ,  3076 ,  3075 ,  3074 ,  3072 ,  3070 ,  3069 ,  3066 ,  3065 ,  3064 ,  3063 ,  3062 ,  3061 ,  3060 ,  3059 ,  3056 ,  3055 ,  3053 ,  3052 ,  3047 ,  3044 ,  3043 ,  3031 ,  3030 ,  3028 ,  3024 ,  3023 ,  3021 ,  3017 ,  3015 ,  3014 ,  3009 ,  3007 ,  3003 ,  3002 ,  3001 ,  3000 ,  2997 ,  2993 ,  2991 ,  2990 ,  2988 ,  2987 ,  2986 ,  2985 ,  2984 ,  2981 ,  2980 ,  2979 ,  2978 ,  2973 ,  2970 ,  2967 ,  2965 ,  2964 ,  2962 ,  2961 ,  2958 ,  2952 ,  2950 ,  2949 ,  2948 ,  2946 ,  2944 ,  2942 ,  2941 ,  2939 ,  2937 ,  2933 ,  2932 ,  2931 ,  2928 ,  2926 ,  2925 ,  2924 ,  2921 ,  2920 ,  2918 ,  2914 ,  2912 ,  2910 ,  2908 ,  2907 ,  2904 ,  2901 ,  2900 ,  2897 ,  2896 ,  2895 ,  2892 ,  2891 ,  2890 ,  2889 ,  2887 ,  2885 ,  2882 ,  2881 ,  2878 ,  2876 ,  2874 ,  2870 ,  2869 ,  2868 ,  2866 ,  2864 ,  2861 ,  2859 ,  2858 ,  2852 ,  2850 ,  2848 ,  2847 ,  2846 ,  2843 ,  2841 ,  2835 ,  2830 ,  2827 ,  2824 ,  2823 ,  2822 ,  2821 ,  2820 ,  2819 ,  2815 ,  2805 ,  2804 ,  2803 ,  2801 ,  2800 ,  2798 ,  2796 ,  2794 ,  2793 ,  2792 ,  2791 ,  2790 ,  2789 ,  2788 ,  2785 ,  2784 ,  2781 ,  2778 ,  2777 ,  2775 ,  2774 ,  2769 ,  2767 ,  2764 ,  2763 ,  2761 ,  2759 ,  2758 ,  2756 ,  2753 ,  2750 ,  2749 ,  2743 ,  2742 ,  2740 ,  2739 ,  2736 ,  2735 ,  2734 ,  2732 ,  2730 ,  2729 ,  2725 ,  2724 ,  2721 ,  2719 ,  2716 ,  2715 ,  2713 ,  2712 ,  2710 ,  2709 ,  2708 ,  2706 ,  2705 ,  2704 ,  2703 ,  2702 ,  2700 ,  2698 ,  2697 ,  2696 ,  2693 ,  2692 ,  2691 ,  2690 ,  2684 ,  2682 ,  2681 ,  2676 ,  2674 ,  2673 ,  2671 ,  2670 ,  2661 ,  2658 ,  2656 ,  2655 ,  2652 ,  2650 ,  2648 ,  2642 ,  2639 ,  2637 ,  2636 ,  2635 ,  2634 ,  2632 ,  2629 ,  2628 ,  2626 ,  2625 ,  2624 ,  2623 ,  2620 ,  2616 ,  2614 ,  2613 ,  2611 ,  2610 ,  2609 ,  2607 ,  2605 ,  2603 ,  2601 ,  2600 ,  2599 ,  2596 ,  2595 ,  2594 ,  2593 ,  2592 ,  2591 ,  2587 ,  2586 ,  2585 ,  2584 ,  2582 ,  2579 ,  2578 ,  2576 ,  2575 ,  2574 ,  2571 ,  2569 ,  2566 ,  2564 ,  2561 ,  2559 ,  2558 ,  2557 ,  2556 ,  2555 ,  2554 ,  2551 ,  2548 ,  2547 ,  2544 ,  2542 ,  2539 ,  2538 ,  2536 ,  2533 ,  2532 ,  2530 ,  2528 ,  2522 ,  2521 ,  2520 ,  2519 ,  2518 ,  2516 ,  2515 ,  2513 ,  2512 ,  2509 ,  2508 ,  2507 ,  2506 ,  2503 ,  2502 ,  2500 ,  2499 ,  2491 ,  2487 ,  2486 ,  2482 ,  2480 ,  2476 ,  2475 ,  2474 ,  2471 ,  2470 ,  2469 ,  2468 ,  2467 ,  2466 ,  2464 ,  2458 ,  2457 ,  2456 ,  2455 ,  2454 ,  2453 ,  2451 ,  2450 ,  2447 ,  2444 ,  2442 ,  2441 ,  2440 ,  2433 ,  2425 ,  2423 ,  2420 ,  2419 ,  2418 ,  2416 ,  2414 ,  2413 ,  2412 ,  2411 ,  2410 ,  2408 ,  2407 ,  2406 ,  2402 ,  2401 ,  2400 ,  2398 ,  2396 ,  2395 ,  2394 ,  2389 ,  2386 ,  2384 ,  2383 ,  2381 ,  2380 ,  2378 ,  2377 ,  2376 ,  2375 ,  2374 ,  2372 ,  2370 ,  2369 ,  2365 ,  2364 ,  2363 ,  2361 ,  2360 ,  2359 ,  2358 ,  2353 ,  2349 ,  2347 ,  2345 ,  2344 ,  2343 ,  2342 ,  2340 ,  2336 ,  2335 ,  2333 ,  2331 ,  2329 ,  2328 ,  2325 ,  2322 ,  2320 ,  2312 ,  2310 ,  2308 ,  2307 ,  2304 ,  2303 ,  2299 ,  2298 ,  2294 ,  2293 ,  2290 ,  2289 ,  2286 ,  2285 ,  2282 ,  2280 ,  2279 ,  2275 ,  2274 ,  2273 ,  2271 ,  2270 ,  2269 ,  2268 ,  2264 ,  2263 ,  2260 ,  2257 ,  2256 ,  2255 ,  2254 ,  2249 ,  2248 ,  2245 ,  2243 ,  2236 ,  2230 ,  2227 ,  2226 ,  2225 ,  2224 ,  2220 ,  2218 ,  2217 ,  2216 ,  2215 ,  2212 ,  2210 ,  2209 ,  2208 ,  2207 ,  2206 ,  2204 ,  2203 ,  2202 ,  2197 ,  2196 ,  2195 ,  2194 ,  2193 ,  2191 ,  2189 ,  2188 ,  2186 ,  2185 ,  2184 ,  2183 ,  2182 ,  2179 ,  2178 ,  2172 ,  2170 ,  2169 ,  2168 ,  2166 ,  2164 ,  2162 ,  2158 ,  2157 ,  2156 ,  2155 ,  2150 ,  2149 ,  2148 ,  2142 ,  2141 ,  2138 ,  2137 ,  2136 ,  2134 ,  2131 ,  2129 ,  2127 ,  2125 ,  2124 ,  2123 ,  2121 ,  2120 ,  2119 ,  2117 ,  2116 ,  2114 ,  2113 ,  2107 ,  2105 ,  2104 ,  2102 ,  2101 ,  2100 ,  2099 ,  2092 ,  2090 ,  2089 ,  2084 ,  2082 ,  2081 ,  2078 ,  2077 ,  2076 ,  2075 ,  2073 ,  2071 ,  2070 ,  2069 ,  2068 ,  2067 ,  2064 ,  2061 ,  2055 ,  2054 ,  2053 ,  2050 ,  2048 ,  2047 ,  2046 ,  2045 ,  2043 ,  2038 ,  2037 ,  2031 ,  2029 ,  2026 ,  2025 ,  2023 ,  2022 ,  2019 ,  2017 ,  2009 ,  2006 ,  2004 ,  2003 ,  1998 ,  1997 ,  1995 ,  1994 ,  1992 ,  1991 ,  1990 ,  1985 ,  1983 ,  1980 ,  1977 ,  1972 ,  1971 ,  1970 ,  1968 ,  1967 ,  1965 ,  1960 ,  1959 ,  1958 ,  1957 ,  1956 ,  1955 ,  1953 ,  1951 ,  1949 ,  1948 ,  1947 ,  1945 ,  1944 ,  1940 ,  1939 ,  1938 ,  1935 ,  1933 ,  1930 ,  1928 ,  1927 ,  1922 ,  1919 ,  1918 ,  1915 ,  1912 ,  1908 ,  1907 ,  1904 ,  1903 ,  1902 ,  1899 ,  1897 ,  1894 ,  1892 ,  1891 ,  1888 ,  1886 ,  1885 ,  1883 ,  1880 ,  1875 ,  1872 ,  1870 ,  1856 ,  1854 ,  1853 ,  1852 ,  1849 ,  1846 ,  1844 ,  1840 ,  1835 ,  1833 ,  1832 ,  1830 ,  1828 ,  1827 ,  1826 ,  1825 ,  1824 ,  1822 ,  1819 ,  1818 ,  1817 ,  1815 ,  1814 ,  1812 ,  1809 ,  1808 ,  1807 ,  1804 ,  1803 ,  1802 ,  1800 ,  1799 ,  1795 ,  1793 ,  1791 ,  1785 ,  1783 ,  1782 ,  1781 ,  1780 ,  1779 ,  1778 ,  1775 ,  1774 ,  1772 ,  1769 ,  1767 ,  1766 ,  1765 ,  1760 ,  1758 ,  1757 ,  1753 ,  1746 ,  1745 ,  1743 ,  1741 ,  1738 ,  1736 ,  1734 ,  1730 ,  1729 ,  1728 ,  1727 ,  1725 ,  1724 ,  1723 ,  1718 ,  1716 ,  1711 ,  1706 ,  1702 ,  1701 ,  1699 ,  1695 ,  1694 ,  1685 ,  1684 ,  1680 ,  1679 ,  1677 ,  1675 ,  1671 ,  1669 ,  1666 ,  1665 ,  1662 ,  1661 ,  1658 ,  1657 ,  1656 ,  1655 ,  1654 ,  1652 ,  1649 ,  1648 ,  1647 ,  1645 ,  1643 ,  1642 ,  1641 ,  1640 ,  1638 ,  1637 ,  1636 ,  1634 ,  1633 ,  1632 ,  1631 ,  1630 ,  1626 ,  1624 ,  1620 ,  1619 ,  1618 ,  1617 ,  1614 ,  1613 ,  1612 ,  1605 ,  1603 ,  1602 ,  1601 ,  1599 ,  1598 ,  1597 ,  1596 ,  1595 ,  1594 ,  1591 ,  1590 ,  1586 ,  1584 ,  1579 ,  1578 ,  1575 ,  1574 ,  1569 ,  1567 ,  1564 ,  1563 ,  1561 ,  1559 ,  1558 ,  1557 ,  1556 ,  1550 ,  1549 ,  1547 ,  1546 ,  1543 ,  1541 ,  1540 ,  1539 ,  1538 ,  1537 ,  1536 ,  1534 ,  1529 ,  1527 ,  1526 ,  1525 ,  1521 ,  1520 ,  1519 ,  1514 ,  1513 ,  1511 ,  1510 ,  1508 ,  1503 ,  1501 ,  1499 ,  1498 ,  1494 ,  1493 ,  1492 ,  1491 ,  1490 ,  1488 ,  1487 ,  1485 ,  1483 ,  1480 ,  1479 ,  1478 ,  1477 ,  1476 ,  1472 ,  1464 ,  1463 ,  1462 ,  1459 ,  1458 ,  1457 ,  1452 ,  1451 ,  1450 ,  1449 ,  1448 ,  1447 ,  1446 ,  1445 ,  1444 ,  1443 ,  1441 ,  1439 ,  1437 ,  1436 ,  1435 ,  1434 ,  1430 ,  1429 ,  1425 ,  1424 ,  1423 ,  1420 ,  1419 ,  1418 ,  1416 ,  1414 ,  1413 ,  1408 ,  1400 ,  1398 ,  1397 ,  1395 ,  1394 ,  1393 ,  1390 ,  1388 ,  1386 ,  1384 ,  1383 ,  1382 ,  1378 ,  1377 ,  1376 ,  1372 ,  1370 ,  1368 ,  1367 ,  1362 ,  1359 ,  1358 ,  1357 ,  1356 ,  1351 ,  1350 ,  1348 ,  1347 ,  1344 ,  1343 ,  1342 ,  1338 ,  1337 ,  1335 ,  1333 ,  1332 ,  1329 ,  1325 ,  1324 ,  1319 ,  1318 ,  1317 ,  1310 ,  1309 ,  1302 ,  1299 ,  1297 ,  1290 ,  1289 ,  1288 ,  1287 ,  1284 ,  1283 ,  1282 ,  1281 ,  1280 ,  1279 ,  1277 ,  1275 ,  1271 ,  1265 ,  1262 ,  1260 ,  1255 ,  1254 ,  1253 ,  1251 ,  1246 ,  1245 ,  1244 ,  1243 ,  1239 ,  1238 ,  1237 ,  1236 ,  1231 ,  1229 ,  1228 ,  1227 ,  1224 ,  1222 ,  1221 ,  1220 ,  1219 ,  1217 ,  1214 ,  1213 ,  1212 ,  1211 ,  1210 ,  1208 ,  1204 ,  1202 ,  1201 ,  1200 ,  1199 ,  1198 ,  1195 ,  1194 ,  1191 ,  1190 ,  1186 ,  1185 ,  1183 ,  1182 ,  1181 ,  1180 ,  1176 ,  1175 ,  1173 ,  1172 ,  1167 ,  1164 ,  1161 ,  1160 ,  1159 ,  1158 ,  1156 ,  1154 ,  1153 ,  1151 ,  1149 ,  1148 ,  1147 ,  1146 ,  1145 ,  1143 ,  1142 ,  1140 ,  1138 ,  1134 ,  1131 ,  1126 ,  1125 ,  1122 ,  1119 ,  1116 ,  1115 ,  1113 ,  1110 ,  1105 ,  1102 ,  1099 ,  1097 ,  1094 ,  1090 ,  1089 ,  1088 ,  1087 ,  1083 ,  1079 ,  1077 ,  1076 ,  1075 ,  1073 ,  1069 ,  1068 ,  1067 ,  1063 ,  1060 ,  1059 ,  1056 ,  1055 ,  1053 ,  1052 ,  1049 ,  1047 ,  1046 ,  1045 ,  1041 ,  1040 ,  1037 ,  1035 ,  1034 ,  1033 ,  1032 ,  1030 ,  1029 ,  1027 ,  1025 ,  1024 ,  1022 ,  1017 ,  1015 ,  1011 ,  1010 ,  1004 ,  1003 ,  1001 ,  1000 ,  999 ,  998 ,  997 ,  996 ,  994 ,  987 ,  978 ,  977 ,  974 ,  972 ,  971 ,  968 ,  967 ,  966 ,  964 ,  963 ,  959 ,  958 ,  956 ,  955 ,  954 ,  952 ,  951 ,  945 ,  944 ,  943 ,  940 ,  938 ,  937 ,  936 ,  933 ,  925 ,  924 ,  922 ,  916 ,  915 ,  914 ,  911 ,  910 ,  907 ,  903 ,  898 ,  894 ,  893 ,  892 ,  890 ,  889 ,  886 ,  884 ,  883 ,  882 ,  881 ,  879 ,  877 ,  874 ,  873 ,  872 ,  871 ,  870 ,  868 ,  866 ,  865 ,  864 ,  861 ,  859 ,  858 ,  857 ,  855 ,  853 ,  841 ,  839 ,  837 ,  835 ,  833 ,  832 ,  830 ,  827 ,  825 ,  824 ,  822 ,  821 ,  820 ,  818 ,  816 ,  814 ,  809 ,  807 ,  806 ,  805 ,  804 ,  803 ,  802 ,  800 ,  799 ,  798 ,  797 ,  795 ,  791 ,  790 ,  789 ,  788 ,  785 ,  784 ,  782 ,  779 ,  776 ,  775 ,  774 ,  772 ,  771 ,  770 ,  769 ,  768 ,  767 ,  765 ,  764 ,  763 ,  761 ,  759 ,  758 ,  754 ,  753 ,  751 ,  749 ,  747 ,  746 ,  745 ,  743 ,  742 ,  740 ,  737 ,  735 ,  730 ,  729 ,  727 ,  726 ,  723 ,  720 ,  719 ,  718 ,  716 ,  715 ,  714 ,  713 ,  712 ,  708 ,  707 ,  706 ,  704 ,  703 ,  702 ,  701 ,  700 ,  698 ,  696 ,  694 ,  692 ,  691 ,  690 ,  688 ,  687 ,  686 ,  684 ,  683 ,  681 ,  680 ,  676 ,  675 ,  674 ,  671 ,  668 ,  667 ,  666 ,  665 ,  664 ,  663 ,  661 ,  655 ,  651 ,  648 ,  647 ,  646 ,  645 ,  643 ,  641 ,  640 ,  638 ,  634 ,  633 ,  630 ,  629 ,  628 ,  627 ,  623 ,  622 ,  620 ,  618 ,  615 ,  613 ,  612 ,  609 ,  604 ,  602 ,  594 ,  592 ,  591 ,  588 ,  582 ,  580 ,  576 ,  572 ,  571 ,  568 ,  564 ,  560 ,  559 ,  555 ,  554 ,  550 ,  549 ,  548 ,  547 ,  545 ,  543 ,  538 ,  537 ,  536 ,  535 ,  532 ,  530 ,  528 ,  527 ,  526 ,  523 ,  521 ,  520 ,  515 ,  511 ,  506 ,  505 ,  504 ,  503 ,  502 ,  501 ,  499 ,  498 ,  494 ,  493 ,  492 ,  490 ,  487 ,  486 ,  485 ,  483 ,  482 ,  481 ,  480 ,  478 ,  477 ,  473 ,  472 ,  470 ,  469 ,  466 ,  465 ,  464 ,  461 ,  460 ,  458 ,  456 ,  454 ,  452 ,  448 ,  447 ,  446 ,  445 ,  444 ,  443 ,  438 ,  437 ,  436 ,  429 ,  428 ,  425 ,  424 ,  420 ,  419 ,  416 ,  409 ,  408 ,  407 ,  406 ,  402 ,  401 ,  400 ,  398 ,  396 ,  395 ,  394 ,  393 ,  392 ,  387 ,  385 ,  384 ,  383 ,  382 ,  381 ,  379 ,  378 ,  377 ,  376 ,  375 ,  365 ,  363 ,  362 ,  360 ,  359 ,  358 ,  353 ,  352 ,  351 ,  350 ,  349 ,  343 ,  341 ,  340 ,  337 ,  335 ,  330 ,  329 ,  328 ,  327 ,  325 ,  324 ,  323 ,  322 ,  320 ,  317 ,  315 ,  314 ,  312 ,  310 ,  307 ,  306 ,  300 ,  299 ,  298 ,  297 ,  292 ,  291 ,  289 ,  288 ,  286 ,  284 ,  283 ,  280 ,  279 ,  278 ,  276 ,  273 ,  271 ,  267 ,  266 ,  265 ,  263 ,  262 ,  261 ,  260 ,  258 ,  257 ,  251 ,  250 ,  248 ,  245 ,  244 ,  243 ,  241 ,  240 ,  237 ,  236 ,  233 ,  231 ,  228 ,  225 ,  221 ,  215 ,  213 ,  212 ,  211 ,  210 ,  208 ,  206 ,  205 ,  202 ,  200 ,  195 ,  191 ,  190 ,  185 ,  183 ,  182 ,  181 ,  179 ,  178 ,  177 ,  175 ,  174 ,  169 ,  165 ,  163 ,  162 ,  161 ,  160 ,  158 ,  155 ,  151 ,  150 ,  149 ,  148 ,  143 ,  142 ,  140 ,  135 ,  133 ,  131 ,  130 ,  127 ,  126 ,  125 ,  124 ,  122 ,  121 ,  120 ,  119 ,  118 ,  117 ,  113 ,  112 ,  111 ,  107 ,  105 ,  103 ,  102 ,  101 ,  100 ,  96 ,  95 ,  93 ,  92 ,  90 ,  87 ,  86 ,  82 ,  81 ,  80 ,  76 ,  74 ,  71 ,  68 ,  66 ,  64 ,  62 ,  59 ,  57 ,  55 ,  50 ,  48 ,  47 ,  44 ,  42 ,  41 ,  40 ,  37 ,  35 ,  31 ,  29 ,  28 ,  25 ,  22 ,  21 ,  20 ,  19 ,  18 ,  16 ,  15 ,  13 ,  12 ,  11 ,  8 ,  7 , 0]

add=[0]*no_clocks
S = s+add
Final_System = [0]*(no_clocks-state_size-max(Alpha))

#Generate equations for update using output bits
# S[15] = S[0]+S[7]
for i in range(no_clocks):
	S[(i+state_size)] =(S[(i+state_size)-20]+S[(i+state_size)-9]+S[(i+state_size)-5]+S[(i+state_size)-3])

for j in range(len(Final_System)):
    for i in Alpha:
    #s1s2z=s1s2s9
    #Here we are going to build two systems; one containing z, the other not
           #Left[j]  = Left[j]+(S[(i+state_size)-(state_size-3)]*S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-6 )]+ S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-10)]*S[(i+state_size)-(state_size-6)])
        Final_System [j] = Final_System[j]+S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-13 )]+ S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-13 )]+ S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-13)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-13)]+(S[(i+state_size+j)-(state_size-13)]+S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-9)]+S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-9)]+S[(i+state_size+j)-(state_size-9)]*S[(i+state_size+j)-(state_size-13)])*O[i+j]
# for i in range(len(Final_System)):
#     for j in range(len(O)):
#         m = eval("x" + str(state_size+j))
#         n = eval('O['+str(j)+']')
#         Final_System[i]=Final_System[i].subs({m:n})
#Define linear substitution variables
stopprecomp = time.perf_counter()
print("precomp",stopprecomp-startprecomp)
start = time.perf_counter()
U = [0]*no_monomials
for i in range(no_monomials):
    U[i] = eval("x" + str(state_size+i))

monomials =[1]*1
for i in range(1,len(U)):
    C = Combinations(s,i)
    no = C.cardinality()
    monomials_temp = [1]*no
    for i in range(no):
        for x in C[i]:
            monomials_temp[i]=monomials_temp[i]*x
    monomials = monomials+monomials_temp

#initialise the dictionary
d = {}
#fill the dictionary with the substituion variables
for i in range(1,len(U)+1):
    d[monomials[i]] = U[i-1]

Final_System_SubWithoutY = [0]*len(Final_System)
Final_System_SubWithY = [0]*len(Final_System)
Final_System_Sub = [0]*len(Final_System)

#For each equation in the system
for i in range(len(Final_System)):

    #NOTE:if you sub for the output bits first, you then only need to do the first line of code as there will no longer be output variables
    #make substitution for monomials without output bits
    Final_System_Sub[i] = sum(Final_System[i].monomial_coefficient(m) * v for m,v in d.items())

F = GF(2)
M = matrix(F,len(Final_System_Sub),len(U))
#Fill coefficient matrix with 1s in posisitons corresponding to initial
# state bits and odd keystream bits
for j in range(len(Final_System_Sub)):
    for i in range(len(U)):
        M[j,i] = Final_System_Sub[j].monomial_coefficient(U[i])
w= vector([0]*len(Final_System_Sub))
# w = [0]*len(Final_System)

# M = transpose(M)

for i in range(len(Final_System)):
#     if (EQ_1[i].monomial_coefficient(1)==1):
    if (Final_System[i].monomial_coefficient(1)==1):      
        w[i]=1
# print(w.parent())
print(M.solve_right(w))
A = M.right_kernel()
# # print(A.dimension())
print(A)
for i in range(A.dimension()):
    print(A[i])



# original_output = sys.stdout
# fp = open("Better_n_15.txt","w")
# sys.stdout = fp

# for i in range(A.dimension()):
#     print(A[i])

stop = time.perf_counter()
print("online",stop-start)


