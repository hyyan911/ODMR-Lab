﻿using CodeHelper;
using HardWares;
using HardWares.相机_CCD_;
using HardWares.端口基类;
using System;
using System.Collections.Generic;
using System.Data.SqlTypes;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace ODMR_Lab
{
    public abstract class DeviceInfoBase<T>
    {

        public bool IsWriting { get; private set; } = false;

        public DeviceConnectInfo ConnectInfo { get; protected set; } = null;

        /// <summary>
        /// 设备名
        /// </summary>
        public T Device { get; set; }

        public abstract DeviceInfoBase<T> CreateDeviceInfo(T device, DeviceConnectInfo connectinfo);

        /// <summary>
        /// 设备关闭事件
        /// </summary>
        public event Action CloseDeviceEvent = null;

        protected abstract void CloseFileAction(FileObject obj);

        /// <summary>
        /// 关闭设备，保存参数，返回参数文件路径
        /// </summary>
        public string CloseDeviceInfoAndSaveParams(out bool IsSucceedClosed)
        {
            if (!IsWriting)
            {
                PortObject dev = Device as PortObject;
                CloseDeviceEvent?.Invoke();
                FileObject obj = new FileObject();
                if (Device is CameraBase)
                {
                    CameraBase c = Device as CameraBase;
                    obj = c.GenerateParamsFile(c);
                }
                dev.Dispose();
                CloseFileAction(obj);
                obj.Descriptions.Add("FileType", "DeviceParamsFile");
                obj.Descriptions.Add("PortType", Enum.GetName(dev.PortType.GetType(), dev.PortType));
                obj.Descriptions.Add("DeviceType", dev.GetType().ToString());
                obj.Descriptions.Add("DeviceInfoType", GetType().ToString());
                if (dev.PortType == PortType.USB)
                {
                    obj.Descriptions.Add("USBConnectName", ConnectInfo.USBName);
                }
                if (dev.PortType == PortType.COM)
                {
                    obj.Descriptions.Add("COMName", ConnectInfo.COMName);
                    obj.Descriptions.Add("BaudRate", ConnectInfo.COMBaud);
                }
                if (!Directory.Exists(Path.Combine(Environment.CurrentDirectory, "DevParamDir")))
                {
                    Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, "DevParamDir"));
                }

                StringBuilder rBuilder = new StringBuilder(Path.Combine(dev.ProductIdentifier + "_" + dev.ProductName));
                foreach (char rinvalidChar in Path.GetInvalidFileNameChars())
                    rBuilder.Replace(rinvalidChar.ToString(), string.Empty);

                string path = Path.Combine(Path.Combine(Environment.CurrentDirectory, "DevParamDir", rBuilder.ToString()));
                obj.SaveToFile(path);
                IsSucceedClosed = true;
                return path + ".userdat";
            }
            else
            {
                MessageLogger.AddLogger("设备", "此设备正在被使用，无法关闭", MessageTypes.Warning, true, true);
                IsSucceedClosed = false;
                return "";
            }
        }

        /// <summary>
        /// 尝试打开所有设备，导入参数
        /// </summary>
        public static List<DeviceInfoBase<T>> ConnectAllAndLoadParams(string ParamFolder, out List<string> unconnectedDeviceinfo)
        {
            if (!Directory.Exists(ParamFolder))
            {
                Directory.CreateDirectory(ParamFolder);
            }
            List<DeviceInfoBase<T>> ConnectDevices = new List<DeviceInfoBase<T>>();
            List<string> result = new List<string>();
            DirectoryInfo info = new DirectoryInfo(ParamFolder);
            FileInfo[] files = info.GetFiles();
            foreach (var item in files)
            {
                if (item.Name.Contains(".userdat"))
                {
                    try
                    {
                        FileObject obj = FileObject.ReadFromFile(item.FullName);

                        if (obj.Descriptions.Values.Contains("DeviceParamsFile"))
                        {
                            //文件是设备参数的文件类型
                            Dictionary<string, Type> devicetypes = PortObject.GetDeviceNamesWithSameBase(typeof(T));
                            foreach (var device in devicetypes)
                            {
                                if (device.Value.ToString() == obj.Descriptions["DeviceType"])
                                {
                                    //找到确定的设备类型
                                    //打开设备
                                    PortType port = (PortType)Enum.Parse(typeof(PortType), obj.Descriptions["PortType"]);
                                    var deviceinstance = Activator.CreateInstance(device.Value);
                                    bool connectresult = false;

                                    DeviceConnectInfo connectinfo = null;

                                    if (port == PortType.USB)
                                    {
                                        connectresult = (deviceinstance as WinUSBOuterInterface).ConnectUSB(obj.Descriptions["USBConnectName"], out Exception exc);
                                        if (connectresult == false)
                                        {
                                            result.Add(((deviceinstance as PortObject).ProductIdentifier + ":\t\t" + obj.Descriptions["USBConnectName"]));
                                            continue;
                                        }
                                        connectinfo = new DeviceConnectInfo(PortType.USB, obj.Descriptions["USBConnectName"]);
                                    }
                                    if (port == PortType.COM)
                                    {
                                        connectresult = (deviceinstance as COMOuterInterface).ConnectCOM(obj.Descriptions["COMName"], int.Parse(obj.Descriptions["BaudRate"]), out Exception exc);
                                        if (connectresult == false)
                                        {
                                            result.Add(((deviceinstance as PortObject).ProductIdentifier + ":\t\t" + obj.Descriptions["COMName"]));
                                            continue;
                                        }
                                        connectinfo = new DeviceConnectInfo(PortType.COM, obj.Descriptions["COMName"], obj.Descriptions["BaudRate"]);
                                    }
                                    if (port == PortType.Socket)
                                    {

                                    }
                                    Type infoType = Type.GetType(obj.Descriptions["DeviceInfoType"]);
                                    DeviceInfoBase<T> devinfo = null;
                                    MainWindow.Handle.Dispatcher.Invoke(() =>
                                    {
                                        devinfo = (Activator.CreateInstance(infoType) as DeviceInfoBase<T>).CreateDeviceInfo((T)deviceinstance, connectinfo);
                                    });
                                    devinfo.Device = (T)deviceinstance;
                                    devinfo.AutoConnectedAction(obj);
                                    ConnectDevices.Add(devinfo);
                                    if (obj.DataCount() != 0)
                                    {
                                        (devinfo.Device as PortObject).LoadParams(obj);
                                    }
                                }
                            }
                        }
                    }
                    catch (Exception) { }
                }
            }
            unconnectedDeviceinfo = result;
            return ConnectDevices;
        }

        protected abstract void AutoConnectedAction(FileObject file);

        /// <summary>
        /// 按照配置文件打开指定设备并导入参数
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns>
        public static DeviceInfoBase<T> OpenAndLoadParams(string filepath, Type deviceType)
        {
            FileObject obj = FileObject.ReadFromFile(filepath);
            //找到确定的设备类型
            //打开设备
            PortType port = (PortType)Enum.Parse(typeof(PortType), obj.Descriptions["PortType"]);
            var deviceinstance = Activator.CreateInstance(deviceType);
            bool connectresult = false;

            DeviceConnectInfo connectinfo = null;

            if (port == PortType.USB)
            {
                connectresult = (deviceinstance as WinUSBOuterInterface).ConnectUSB(obj.Descriptions["USBConnectName"], out Exception exc);
                if (connectresult == false)
                {
                    throw new Exception("设备连接失败");
                }
                connectinfo = new DeviceConnectInfo(PortType.USB, obj.Descriptions["USBConnectName"]);
            }
            if (port == PortType.COM)
            {
                connectresult = (deviceinstance as COMOuterInterface).ConnectCOM(obj.Descriptions["COMName"], int.Parse(obj.Descriptions["BaudRate"]), out Exception exc);
                if (connectresult == false)
                {
                    throw new Exception("设备连接失败");
                }
                connectinfo = new DeviceConnectInfo(PortType.COM, obj.Descriptions["COMName"], obj.Descriptions["BaudRate"]);
            }
            if (port == PortType.Socket)
            {

            }

            Type infoType = Type.GetType(obj.Descriptions["DeviceInfoType"]);
            DeviceInfoBase<T> devinfo = null;
            MainWindow.Handle.Dispatcher.Invoke(() =>
            {
                devinfo = (Activator.CreateInstance(infoType) as DeviceInfoBase<T>).CreateDeviceInfo((T)deviceinstance, connectinfo);
            });
            devinfo.Device = (T)deviceinstance;
            devinfo.AutoConnectedAction(obj);
            if (obj.DataCount() != 0)
            {
                (devinfo.Device as PortObject).LoadParams(obj);
            }

            return devinfo;
        }

        /// <summary>
        /// 使用设备,向设备写信息，当设备处于使用状态时其他无关操作会导致报错
        /// </summary>
        public void BeginUse()
        {
            IsWriting = true;
        }

        /// <summary>
        /// 停止向设备写信息，当设备处于使用状态时其他无关操作会导致报错
        /// </summary>
        public void EndUse()
        {
            IsWriting = false;
        }
    }
}
